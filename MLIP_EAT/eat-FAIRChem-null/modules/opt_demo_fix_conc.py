#!/usr/bin/env python3
import sys
from pathlib import Path
import math
import numpy as np
import torch
from ase import Atoms
from ase.data import atomic_numbers
import logging
from datetime import datetime
import os, random

# FairChem
from fairchem.core import FAIRChemCalculator, pretrained_mlip
from .energy_above_hull import HullEnergyCalculator
from .minimizer import minimizer


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# Random seed setup since FairChem sets the seed and that makes the random initial EAT weights the same each iteration
def save_rng_state():
    st = {
        "py": random.getstate(),
        "np": np.random.get_state(),
        "torch_cpu": torch.get_rng_state(),
    }
    if torch.cuda.is_available():
        st["torch_cuda"] = torch.cuda.get_rng_state_all()
    return st

def load_rng_state(st):
    random.setstate(st["py"])
    np.random.set_state(st["np"])
    torch.set_rng_state(st["torch_cpu"])
    if torch.cuda.is_available() and "torch_cuda" in st:
        torch.cuda.set_rng_state_all(st["torch_cuda"])
rng_state = save_rng_state()




# ============================================================================
# LOGGING SETUP
# ============================================================================

log_dir = Path("logs")
log_dir.mkdir(exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = log_dir / f"optimization_combo_{timestamp}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

log.info("🚀 Combo Optimization started")
log.info(f"Log file: {log_file}")

# ============================================================================
# OPTIMIZATION HISTORY STORAGE (FOR PLOTTING)
# ============================================================================

history = {
    "iter": [],
    "entropy": [],
    "cost": [],
    "total_energy": [],        
    "formation_energy": [],
    "energy_above_hull": [],
    "cell_volume": [],   # Å^3
    "force_norm": [],    # |grad_frac| per-atom
}


anneal_boundaries = []
global_iter = 0

# ============================================================================
# AUTOGRAD-SAFE LATTICE + FRACTIONAL POSITION SETUP
# ============================================================================

def ucparams_to_lattice(ucparams: torch.Tensor) -> torch.Tensor:
    """
    Converts (a,b,c,alpha,beta,gamma) -> lattice matrix H with rows as lattice vectors.
    Angles must be in radians.
    """
    a, b, c, alpha, beta, gamma = ucparams
    zero = torch.zeros((), device=ucparams.device, dtype=ucparams.dtype)

    e1 = torch.stack([a, zero, zero])
    e2 = torch.stack([
        b * torch.cos(gamma),
        b * torch.sin(gamma),
        zero,
    ])

    e3_0 = c * torch.cos(beta)
    e3_1 = c * (torch.cos(alpha) - torch.cos(beta) * torch.cos(gamma)) / (torch.sin(gamma) + 1e-12)
    e3_2 = torch.sqrt(torch.clamp(c**2 - e3_0**2 - e3_1**2, min=1e-12))
    e3 = torch.stack([e3_0, e3_1, e3_2])

    H = torch.stack((e1, e2, e3), dim=0)  # ROWS as lattice vectors
    return H


def cell_to_ucparams(cell):
    """
    cell is 3x3 with rows as ASE cell vectors.
    This returns (a,b,c,alpha,beta,gamma) in radians.
    """
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    alpha = math.acos(np.dot(cell[1], cell[2]) / (b * c))
    beta  = math.acos(np.dot(cell[0], cell[2]) / (a * c))
    gamma = math.acos(np.dot(cell[0], cell[1]) / (a * b))
    return np.array([a, b, c, alpha, beta, gamma])


def jacobian_ucparams_to_lattice(ucparams: torch.Tensor) -> torch.Tensor:
    """
    Returns J with shape (3,3,6) where J[i,j,k] = dH[i,j]/d(ucparams[k]).
    Uses autograd on ucparams_to_lattice.
    """
    H = ucparams_to_lattice(ucparams)
    J = torch.zeros(3, 3, 6, dtype=ucparams.dtype, device=ucparams.device)

    for i in range(3):
        for j in range(3):
            gij = torch.autograd.grad(H[i, j], ucparams, retain_graph=True)[0]
            J[i, j, :] = gij
    return J

# ============================================================================
# FAIRCHEM EAT HELPERS
# ============================================================================

def attach_eat_fairchem(atoms, eat_weights_np):
    """
    Attach per-atom EAT weights to ASE atoms for FAIRChem.
    eat_weights_np: (n_atoms, 100)
    """
    atoms.new_array("eat_weights", np.asarray(eat_weights_np))
    return atoms


def set_enable_eat_grad_fairchem(atoms, flag: bool):
    """
    FAIRChem hook: AtomicData.from_ase reads this flag.
    """
    atoms.info["enable_eat_grad"] = bool(flag)
    return atoms

# ============================================================================
# FAIRCHEM CALCULATOR
# ============================================================================

predictor = pretrained_mlip.get_predict_unit("uma-s-1")
calc = FAIRChemCalculator(predictor, task_name="omat")

# Load random seed after fairchem defined so EAT weights are random each run
load_rng_state(rng_state)


def get_E_F_S(atoms, eat_weights_np=None, enable_eat_grad=False):
    if eat_weights_np is not None:
        attach_eat_fairchem(atoms, eat_weights_np)
    if enable_eat_grad:
        set_enable_eat_grad_fairchem(atoms, True)

    atoms.calc = calc
    E = atoms.get_potential_energy()
    F = atoms.get_forces()
    S = atoms.get_stress(voigt=False)
    return E, F, S

# ============================================================================
# USER INPUT
# ============================================================================

selected_elements = ["Cu", "Pd"]
x_Cu = 0.5
a = x_Cu * 3.61 + (1-x_Cu) * 3.89
# Lattice vectors are ROWS
lattice_vectors = [
    [a, 0.00, 0.00],
    [0.00, a, 0.00],
    [0.00, 0.00, a],
]
symbols = ["Nb", "Nb", "Nb", "Nb"]
frac_coords_init = [
    [0.0, 0.0, 0.0],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.5],
]
supercell_size = (2, 2, 1)

eta_init = 10.0
n_anneals = 1000
anneal_mult = 1.0
BFGS_hist = 5
maxit_per_eta = 100
tol_pg = 1e-2
tol_cost = 1e-3
nIterations_tol = 5
length_rel_delta = 0.5
angle_deg_margin = 30.0
angle_min_deg = 40.0
angle_max_deg = 140.0

Tsallis_q = 2.0

gamma = 10              # entropy cost sharpness
tau_mult = 0.8          # entropy upper bound decay per anneal

comp_target = torch.tensor([x_Cu, 1-x_Cu], device=device) # target composition of Cu and Pd we want to hit


# ============================================================================
# PREPARE HULL ENERGY CALCULATOR
# ===========================================================================
DB_JSON_PATH = "temp_db.json"  # Adjust path as needed
hullcalc = HullEnergyCalculator(DB_JSON_PATH, selected_elements)

def calc_hull(x):
    composition = torch.sum(x, dim=0)
    composition = composition / torch.sum(composition)
    composition = composition.detach().cpu().numpy()
    E_hull, grad_c = hullcalc.get_energy_and_gradient(composition) # gives energy per atom
    grad_x = torch.tensor(grad_c, dtype=x.dtype, device=x.device).unsqueeze(0).expand_as(x) / torch.sum(x) # divide by total number of atoms (sum over x)
    return E_hull, grad_x.flatten()

# ============================================================================
# BUILD BASE STRUCTURE
# ============================================================================

atoms_base = Atoms(
    symbols=symbols,
    scaled_positions=frac_coords_init,
    cell=lattice_vectors,
    pbc=True
)
atoms_base = atoms_base * supercell_size

n_atoms = len(atoms_base)
symbols_all = atoms_base.get_chemical_symbols()

frac0_np = atoms_base.get_scaled_positions()
cell0_np = atoms_base.cell.array.copy()

uc0_np = cell_to_ucparams(cell0_np)
uc0 = torch.tensor(uc0_np, dtype=torch.float32, device=device)

# ============================================================================
# DUAL ELEMENT BASIS: FAIRCHEM
# ============================================================================

n_fairchem_elements = 100
fairchem_element_indices = torch.tensor(
    [atomic_numbers[elem] for elem in selected_elements],
    dtype=torch.long,
    device=device,
)
n_elements = len(selected_elements)
fairchem_indices = fairchem_element_indices.unsqueeze(0).expand(n_atoms, -1)

# ============================================================================
# SIMPLEX PROJECTION
# ============================================================================

def project_simplex_rows(X):
    """
    Row-wise Euclidean projection onto simplex.
    """
    X_sorted, _ = torch.sort(X, dim=1, descending=True)
    cumsum = torch.cumsum(X_sorted, dim=1)
    rho = torch.arange(1, X.shape[1] + 1, device=X.device).view(1, -1)

    cond = X_sorted - (cumsum - 1) / rho > 0
    k = torch.clamp(cond.sum(dim=1), min=1)
    idx = (k - 1).long()

    theta = (cumsum[torch.arange(X.shape[0]), idx] - 1) / k
    X_proj = torch.clamp(X - theta.view(-1, 1), min=0.0)
    return X_proj


def project_flat_all(z_flat, N, S, uc_lower, uc_upper):
    """
    Projects:
      - EAT rows -> simplex
      - ucparams -> clamp to bounds
      - fractional positions -> clamp into (0,1)
    """
    n_x = N * S

    x_flat = z_flat[:n_x]
    uc = z_flat[n_x:n_x + 6]
    frac = z_flat[n_x + 6:].view(N, 3)

    X = x_flat.view(N, S)
    X_proj = project_simplex_rows(X).view(-1)

    uc_proj = torch.clamp(uc, uc_lower, uc_upper)

    return torch.cat([X_proj, uc_proj, frac.view(-1)], dim=0)

# ============================================================================
# ENTROPY + GRAD
# ============================================================================

def entropy(X, Tsallis_q, eps=1e-10):
    """
    Returns:
      ent: scalar
      grad_ent: same shape as X
    """
    if Tsallis_q == 1.0:
        ent = -(X * torch.log(X + eps)).sum()
        grad_ent = -(torch.log(X + eps) + 1.0)
    else:
        ent = torch.sum(1.0 - torch.sum(X**Tsallis_q, dim=-1)) / (Tsallis_q - 1.0)
        grad_ent = -Tsallis_q * X**(Tsallis_q - 1.0) / (Tsallis_q - 1.0)

    return ent, grad_ent

# ============================================================================
# FORMATION ENERGY COST
# ============================================================================

def softplus(E, beta):
    """
    Softplus formation energy cost function and gradient factor.
    """
    f = math.log(1 + math.exp(beta * E)) / beta
    g = math.exp(beta * E) / (1 + math.exp(beta * E))
    return f, g

# ============================================================================
# FORMATION ENERGY SETUP (JSON REFERENCE — CORRECT)
# ============================================================================

from formation_energy_helpers import FormationEnergyCalculator
import json

# 🔽 ABSOLUTE PATH TO YOUR JSON FILE 🔽
REFERENCE_JSON_PATH = Path("reference_energies.json")

if not REFERENCE_JSON_PATH.exists():
    raise RuntimeError(f"Reference energy JSON NOT FOUND: {REFERENCE_JSON_PATH.resolve()}")

with open(REFERENCE_JSON_PATH, "r") as f:
    ref_energies = json.load(f)

# ✅ Build formation-energy calculator from JSON ONLY
formation_calc = FormationEnergyCalculator(
    reference_energies_path=REFERENCE_JSON_PATH
)

# ✅ Build tensor in *design-basis order* for analytic gradients
ref_energies_tensor = torch.tensor(
    [ref_energies[e] for e in selected_elements],
    dtype=torch.float32,
    device=device,
)

# ============================================================================
# SHARED METRICS
# ============================================================================

last_metrics = {
    "entropy_per_atom": None,
    "entropy_total": None,
    "cost": None,
    "total_energy": None,
    "formation_energy": None,
    "energy_above_hull": None,       
    "cell_volume": None,
    "force_norm": None,
}

# --- NEW: Early-stop machinery for PBFGS plateau convergence ---

class PBFGSPlateauReached(Exception):
    """Raised when the PBFGS cost has plateaued for nIterations_tol steps."""
    pass

# These will be reset at each annealing stage
plateau_count = 0
last_cost_for_plateau = None
z_last = None   # last iterate seen by the callback

class PBFGSForceConverged(Exception):
    """Raised when PBFGS max force falls below a specified tolerance."""
    pass

# Global force tolerance for PBFGS; None = disabled
FORCE_TOL = None

# Whether to use cost-plateau early stopping in pbfgs_callback
USE_COST_PLATEAU = True


# ============================================================================
# OBJECTIVE: ENTROPY + FAIRCHEM **FORMATION** ENERGY
# Layout: z = [x_raw, uc_raw, frac_raw]
# ============================================================================

def make_fun_jac(eta, tau, relax):

    n_x = n_atoms * n_elements
    n_uc = 6

    def fun_jac(z_flat):

        # Unpack with fresh leaf tensors so grads are clean
        x_raw = z_flat[:n_x].detach().clone().requires_grad_(True)
        uc    = z_flat[n_x:n_x + n_uc].detach().clone().requires_grad_(True)
        frac  = z_flat[n_x + n_uc:].detach().clone().view(n_atoms, 3).requires_grad_(True)

        # -------------------------
        # EAT weights (project → X)
        # -------------------------
        X = project_simplex_rows(x_raw.view(n_atoms, n_elements))
        entropy_total, grad_entropy = entropy(X, Tsallis_q)
        grad_entropy_flat = grad_entropy.reshape(-1)
        entropy_cost, grad_entropy_cost = softplus(entropy_total - n_atoms*tau, gamma)
        grad_entropy_cost = grad_entropy_cost * grad_entropy_flat

        # Build FairChem EAT (n_atoms, 100)
        eat_fairchem = torch.zeros(
            (n_atoms, n_fairchem_elements), device=device, dtype=X.dtype
        )
        eat_fairchem.scatter_add_(1, fairchem_indices, X)

        # -------------------------
        # Lattice
        # -------------------------
        cell = ucparams_to_lattice(uc)

        # -------------------------
        # ===== FAIRCHEM ENERGY + MANUAL GRADS =====
        # -------------------------
        atoms = atoms_base.copy()

        atoms.set_scaled_positions(frac.detach().cpu().numpy())
        atoms.set_cell(cell.detach().cpu().numpy(), scale_atoms=True)

        eat_np = eat_fairchem.detach().cpu().numpy()
        E, F_cart, S = get_E_F_S(atoms, eat_weights_np=eat_np, enable_eat_grad=True)

        if "eat_grad" not in calc.results:
            raise RuntimeError(
                "FAIRChem did not return 'eat_grad'. "
                "Your fairchem build must support enable_eat_grad."
            )

        # full EAT gradient in 100-element basis
        eat_grad_full = torch.tensor(
            calc.results["eat_grad"],
            device=device,
            dtype=X.dtype
        )  # (n_atoms, 100)

        # map to selected-elements gradient: gather columns at fairchem_indices
        grad_X_fairchem = eat_grad_full.gather(1, fairchem_indices)  # (n_atoms, n_elements)
        grad_x_fairchem = grad_X_fairchem.reshape(-1)

        # fractional gradient from forces
        F_t = torch.tensor(F_cart, device=device, dtype=X.dtype)
        grad_frac_fairchem = -(F_t @ cell.T).reshape(-1)

        # stress -> dE/dH
        V = atoms.get_volume()
        S_t = torch.tensor(S, device=device, dtype=X.dtype)
        HinvT = torch.inverse(cell).T
        dE_dH = V * (HinvT @ S_t)

        # chain rule to ucparams
        J = jacobian_ucparams_to_lattice(uc)  # (3,3,6)
        grad_uc_fairchem = torch.zeros(6, device=device, dtype=X.dtype)
        for k in range(6):
            grad_uc_fairchem[k] = torch.sum(dE_dH * J[:, :, k])

        E_t = torch.tensor(E, device=device, dtype=X.dtype)

        # -------------------------
        # ===== FORMATION ENERGY (VALUE + GRADS) =====
        # -------------------------
        # composition fractions from design-basis X
        comp_frac = X.mean(dim=0)  # (n_elements,)

        composition = {
            selected_elements[i]: float(comp_frac[i].item())
            for i in range(n_elements)
            if comp_frac[i].item() > 1e-8
        }

        # formation energy per atom
        E_form = formation_calc.compute_formation_energy(
            total_energy=E_t,
            composition=composition,
            num_atoms=n_atoms,
        )

        # Energy of hull
        E_hull, grad_x_hull = calc_hull(X)

        # Energy above hull
        E_above_hull = E_form - E_hull

        # d(E/N)/d*
        invN = 1.0 / float(n_atoms)
        grad_x_energy_part = invN * grad_x_fairchem
        grad_uc_energy_part = invN * grad_uc_fairchem
        grad_frac_energy_part = invN * grad_frac_fairchem

        # d(- sum f_j Eref_j)/dX_aj = -(1/N) * Eref_j
        grad_X_comp = (-invN) * ref_energies_tensor.to(dtype=X.dtype).view(1, -1).expand(n_atoms, -1)
        grad_x_comp = grad_X_comp.reshape(-1)

        # total formation-energy gradients
        grad_x_form = grad_x_energy_part + grad_x_comp
        grad_uc_form = grad_uc_energy_part
        grad_frac_form = grad_frac_energy_part

        # total energy above hull gradients
        grad_x_above_hull = grad_x_form - grad_x_hull
        grad_uc_above_hull = grad_uc_form
        grad_frac_above_hull = grad_frac_form

        # composition penalty
        comp_penalty = 5*n_atoms/2 * torch.sum((comp_frac - comp_target)**2)
        comp_grad = torch.zeros_like(X)
        for i in range(len(selected_elements)):
            comp_grad[:, i] = 5*(comp_frac[i] - comp_target[i])
        comp_grad = comp_grad.flatten()

        # -------------------------
        # COMBO COST (use formation energy)
        # -------------------------
        cost = eta * entropy_cost + E_above_hull + comp_penalty
        grad_x =  eta * grad_entropy_cost + grad_x_above_hull + comp_grad
        grad_uc = grad_uc_form
        grad_f  = grad_frac_form
        if relax:
            grad_x = grad_x * 0.0       # no EAT optimization during relaxation
        else:
            grad_uc = grad_uc * 0.0
            grad_f = grad_f * 0.0

        grad = torch.cat([grad_x, grad_uc, grad_f], dim=0)

        if torch.isnan(grad).any() or torch.isinf(grad).any():
            raise RuntimeError("NaN/Inf detected in gradient.")

        with torch.no_grad():
            cell_volume = float(atoms.get_volume())

            # Fractional gradient norm per atom (for diagnostics)
            gf = grad_f.view(n_atoms, 3)
            force_norm = float(torch.norm(gf, dim=1).mean().item())

            # NEW: true max Cartesian force (for PBFGS force-based stopping)
            F_t = torch.tensor(F_cart, device=device, dtype=X.dtype)
            fmax = float(torch.norm(F_t, dim=1).max().item())

            last_metrics["cell_volume"] = cell_volume
            last_metrics["force_norm"] = force_norm
            last_metrics["fmax"] = fmax
            last_metrics["entropy_total"] = entropy_total.detach().item()
            last_metrics["entropy_per_atom"] = (entropy_total / n_atoms).detach().item()
            last_metrics["cost"] = cost.detach().item()
            last_metrics["total_energy"] = float(E)  # total energy
            last_metrics["formation_energy"] = E_form.detach().item()
            last_metrics["energy_above_hull"] = float(E_above_hull.detach().item())

        return cost.detach(), grad.detach()

    return fun_jac

# ============================================================================
# CALLBACK
# ============================================================================

def pbfgs_callback(z, f, err, k, converged, *args):
    global global_iter
    global plateau_count, last_cost_for_plateau, z_last, FORCE_TOL, USE_COST_PLATEAU

    # store latest iterate in case we early-stop
    z_last = z.detach().clone()

    f_val = f.item() if torch.is_tensor(f) else float(f)

    ent = last_metrics.get("entropy_per_atom", None)
    cost = last_metrics.get("cost", None)
    total_energy = last_metrics.get("total_energy", None)
    formation_energy = last_metrics.get("formation_energy", None)
    energy_above_hull = last_metrics.get("energy_above_hull", None)
    cell_volume = last_metrics.get("cell_volume", None)
    force_norm = last_metrics.get("force_norm", None)
    fmax = last_metrics.get("fmax", None)  # NEW

    history["iter"].append(global_iter)
    history["entropy"].append(ent if ent is not None else np.nan)
    history["cost"].append(cost if cost is not None else np.nan)
    history["total_energy"].append(total_energy if total_energy is not None else np.nan)
    history["formation_energy"].append(formation_energy if formation_energy is not None else np.nan)
    history["energy_above_hull"].append(energy_above_hull if energy_above_hull is not None else np.nan)
    history["cell_volume"].append(cell_volume if cell_volume is not None else np.nan)
    history["force_norm"].append(force_norm if force_norm is not None else np.nan)

    # --- plateau detection on cost (can be disabled) ---
    if USE_COST_PLATEAU and cost is not None:
        if last_cost_for_plateau is not None:
            if abs(cost - last_cost_for_plateau) <= tol_cost:
                plateau_count += 1
            else:
                plateau_count = 0
        last_cost_for_plateau = cost

        if plateau_count >= nIterations_tol:
            log.info(
                f"[PBFGS iter {k:03d}] Plateau convergence reached: "
                f"|Δcost| <= {tol_cost} for {plateau_count} steps. "
                f"Stopping PBFGS early for this anneal."
            )
            raise PBFGSPlateauReached()


    # --- NEW: force-based stopping when FORCE_TOL is set ---
    if FORCE_TOL is not None and fmax is not None:
        if fmax <= FORCE_TOL:
            log.info(
                f"[PBFGS iter {k:03d}] Force convergence reached: "
                f"fmax = {fmax:.3e} <= FORCE_TOL = {FORCE_TOL:.3e}. Stopping PBFGS."
            )
            raise PBFGSForceConverged()

    global_iter += 1

    if ent is None or cost is None:
        log.info(
            f"[PBFGS iter {k:03d}] f={f_val:.6f}, PG-error={err:.3e}"
        )
        return

    log.info(
        f"[PBFGS iter {k:03d}] "
        f"cost={cost:.6f}, E_above_hull={energy_above_hull:.6f} eV/atom, "
        f"V={cell_volume:.2f} Å^3, fmax={fmax:.3e}, "
        f"ent/atom={ent:.3e}, PG-error={err:.3e}"
    )




# ============================================================================
# PBFGS DRIVER
# ============================================================================

# bounds on ucparams
uc_lower = uc0.clone()
uc_upper = uc0.clone()

for i in range(3):
    uc_lower[i] = uc0[i] * (1.0 - length_rel_delta)
    uc_upper[i] = uc0[i] * (1.0 + length_rel_delta)

angle_margin_rad = math.radians(angle_deg_margin)
angle_min_rad = math.radians(angle_min_deg)
angle_max_rad = math.radians(angle_max_deg)

for i in range(3, 6):
    uc_lower[i] = max(uc0[i].item() - angle_margin_rad, angle_min_rad)
    uc_upper[i] = min(uc0[i].item() + angle_margin_rad, angle_max_rad)

uc_lower = uc_lower.to(device)
uc_upper = uc_upper.to(device)

n_x = n_atoms * n_elements
n_uc = 6

# init EAT in design basis (selected_elements) -- NEED TO MAKE RANDOM EACH RUN - right now each indep run gives same x0!
dist = torch.distributions.Dirichlet(torch.ones(n_elements, device=device))
x0 = dist.sample((n_atoms,))
x0 = project_simplex_rows(x0).view(-1)

# init frac
f0 = torch.tensor(frac0_np, dtype=torch.float32, device=device)

# initial z
z0 = torch.cat([x0, uc0, f0.reshape(-1)], dim=0)
tau_init = 1.0 - (1.0 / n_atoms)
for anneal in range(n_anneals):
    anneal_boundaries.append(global_iter)

    eta   = eta_init   * (anneal_mult ** anneal)
    tau = tau_init * (tau_mult ** anneal)

    log.info("")
    log.info("=" * 80)
    log.info(
        f"ANNEALING STEP {anneal + 1}/{n_anneals} | "
        f"eta = {eta:.6e} | tau = {tau:.6e}"
    )
    log.info("=" * 80)

    fun_jac = make_fun_jac(eta, tau, relax=False)

    # --- NEW: reset plateau tracking for this anneal ---
    plateau_count = 0
    last_cost_for_plateau = None

    opt = minimizer(
        fun_jac=fun_jac,
        x0=z0,
        args=(),
        project=project_flat_all,
        proj_args=(n_atoms, n_elements, uc_lower, uc_upper),
        callback=pbfgs_callback,
        BFGS_hist=BFGS_hist,
        maxit=maxit_per_eta,
        tol=tol_pg,
        max_line_search=5,
    )

    try:
        obj, z0 = opt.PBFGS()
        X = z0[:n_x].view(n_atoms, n_elements)
        ent = entropy(X, Tsallis_q)[0]
        if ent < 1e-4:
            X_rounded = torch.zeros_like(X)
            max_indices = torch.argmax(X, dim=1)
            X_rounded.scatter_(1, max_indices.unsqueeze(1), 1.0)
            X_rounded = X_rounded.view(-1)
            z0 = torch.cat([X_rounded, z0[n_x:]], dim=0)
            log.info("✅ Entropy sufficiently small. Rounding weights and moving onto final relaxation.")
            break
    except PBFGSPlateauReached:
        log.info("PBFGS plateau convergence criterion reached; moving to next anneal.")
        # use the last iterate seen by the callback as the new starting point
        if z_last is not None:
            z0 = z_last.detach().clone()
        else:
            log.warning("z_last is None after plateau; keeping previous z0.")
        X = z0[:n_x].view(n_atoms, n_elements)
        ent = entropy(X, Tsallis_q)[0]
        if ent < 1e-4:
            X_rounded = torch.zeros_like(X)
            max_indices = torch.argmax(X, dim=1)
            X_rounded.scatter_(1, max_indices.unsqueeze(1), 1.0)
            print(X_rounded)
            X_rounded = X_rounded.view(-1)
            z0 = torch.cat([X_rounded, z0[n_x:]], dim=0)
            log.info("✅ Entropy sufficiently small. Rounding weights and moving onto final relaxation.")
            break




# ============================================================================
# FINAL PBFGS RELAXATION (α = 0, allow EAT, uc, frac to move)
# ============================================================================

eta_final = 0.0 # syntropization not needed - assumed material is already syntropized
tau_relax = 0.0

log.info("")
log.info("=" * 80)
log.info("🔧 FINAL PBFGS RELAXATION (α = 0, force-based stop)")
log.info("=" * 80)

fun_jac_relax = make_fun_jac(eta_final, tau_relax, relax=True)

# --- CRITICAL: final relax uses ONLY force stopping ---
maxit_final = 500
FORCE_TOL = 1e-2       # your fmax threshold
USE_COST_PLATEAU = False # disable cost plateau for this run

plateau_count = 0
last_cost_for_plateau = None

opt_relax = minimizer(
    fun_jac=fun_jac_relax,
    x0=z0,
    args=(),
    project=project_flat_all,
    proj_args=(n_atoms, n_elements, uc_lower, uc_upper),
    callback=pbfgs_callback,
    BFGS_hist=BFGS_hist,
    maxit=maxit_final,  # still a hard cap
    tol=0.0,              # effectively disables PG-based convergence
    max_line_search=5,
)

try:
    obj_relax, z0 = opt_relax.PBFGS()
except PBFGSForceConverged:
    log.info("✅ Final PBFGS relaxation: force convergence reached.")
    if z_last is not None:
        z0 = z_last.detach().clone()
except PBFGSPlateauReached:
    # Should not happen with USE_COST_PLATEAU=False, but kept as guard
    log.info("⚠️ PlateauReached during final relaxation (unexpected).")
    if z_last is not None:
        z0 = z_last.detach().clone()

FORCE_TOL = None
USE_COST_PLATEAU = True   # reset for any future runs






# ============================================================================
# FINAL RESULT LOGGING
# ============================================================================

log.info("")
log.info("✅ OPTIMIZATION COMPLETE")
log.info("=" * 80)

x_final = z0[:n_x]
uc_final = z0[n_x:n_x + n_uc]
frac_final = z0[n_x + n_uc:].view(n_atoms, 3)

cell_final = ucparams_to_lattice(uc_final).detach().cpu().numpy()
frac_final_np = frac_final.detach().cpu().numpy()

# Rebuild X, Tc EAT, and FairChem EAT at optimum
Xf = project_simplex_rows(x_final.view(n_atoms, n_elements))

eat_fairchem_final = torch.zeros((n_atoms, n_fairchem_elements), device=device)
eat_fairchem_final.scatter_add_(1, fairchem_indices, Xf)

# Final FairChem energy evaluation
atoms_fin = atoms_base.copy()
atoms_fin.set_scaled_positions(frac_final_np)
atoms_fin.set_cell(cell_final, scale_atoms=True)
E_fin, _, _ = get_E_F_S(
    atoms_fin,
    eat_weights_np=eat_fairchem_final.detach().cpu().numpy(),
    enable_eat_grad=False
)

# Final formation energy per atom
with torch.no_grad():
    comp_frac_final = Xf.mean(dim=0)
    composition_final = {
        selected_elements[i]: float(comp_frac_final[i].item())
        for i in range(n_elements)
        if comp_frac_final[i].item() > 1e-8
    }
    E_fin_t = torch.tensor(E_fin, device=device, dtype=Xf.dtype)
    E_form_fin = formation_calc.compute_formation_energy(
        total_energy=E_fin_t,
        composition=composition_final,
        num_atoms=n_atoms,
    ).item()
    E_hull_fin, _ = calc_hull(Xf)
    E_above_hull_fin = E_form_fin - E_hull_fin

log.info(f"FINAL E_total  = {E_fin:.8f} eV")
log.info(f"FINAL E_above_hull   = {E_above_hull_fin:.8f} eV/atom")

log.info("")
log.info("FINAL UNIT CELL (3x3):")
for row in cell_final:
    log.info("  " + " ".join(f"{x: .8f}" for x in row))

log.info("")
log.info("FINAL FRACTIONAL POSITIONS:")
for i, xyz in enumerate(frac_final_np):
    log.info(f"Atom {i:03d}:  {xyz[0]: .8f}  {xyz[1]: .8f}  {xyz[2]: .8f}")

# Optional: print Tc-basis EAT weights
log.info("")
log.info("FINAL EAT WEIGHTS (rows=atoms, cols=Tc element types with nonzero weights)")
log.info("=" * 100)

header = "Atom \\ Element | " + " ".join([f"{elem:>6s}" for elem in selected_elements])
log.info(header)
log.info("-" * len(header))

for i in range(n_atoms):
    row_vals = Xf[i].detach().cpu().numpy()
    row_str = f"{i:14d} | " + " ".join([f"{v:6.3f}" for v in row_vals])
    log.info(row_str)

log.info("=" * 100)

log.info("")
log.info("✅ FULL RESULTS WRITTEN TO LOG FILE")
log.info("=" * 80)





# Add new material to hull database
fractions = []
label = ''
for elem, comp in zip(selected_elements, comp_frac_final):
    label += f'{elem}{comp.item():.3f} '
    fractions.append(comp.item())
hullcalc.add_phase_to_db(label, selected_elements, fractions, E_form_fin)






# ============================================================================
# PLOTTING: dashboard of Tc, E_form, cost, entropy, volume, grad + final cell
# 2x4 grid
# ============================================================================

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# optional ASE render helper
try:
    from ase.visualize.plot import plot_atoms
    ASE_PLOT_OK = True
except Exception:
    ASE_PLOT_OK = False

iters = np.array(history["iter"])
ent_vals = np.array(history["entropy"])
cost_vals = np.array(history["cost"])
E_vals = np.array(history["total_energy"])
Eform_vals = np.array(history["formation_energy"])
energy_above_hull_vals = np.array(history["energy_above_hull"])
vol_vals = np.array(history.get("cell_volume", []))
force_norm_vals = np.array(history.get("force_norm", []))

if len(iters) > 0:

    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    ax = axes.flatten()

    def add_anneal_lines(a):
        for b in anneal_boundaries:
            a.axvline(b, linestyle="--", linewidth=0.8)

    # ========== 0) Cost ==========
    ax[0].plot(iters, cost_vals, "k-o", ms=3)
    add_anneal_lines(ax[0])
    ax[0].set_title("Cost (entropy + E_form)")
    ax[0].set_xlabel("Iteration")
    ax[0].set_ylabel("Cost")
    ax[0].grid(True, alpha=0.3)

    # ========== 2) Cell Volume ==========
    if len(vol_vals) == len(iters):
        ax[1].plot(iters, vol_vals, "b-o", ms=3)
        add_anneal_lines(ax[1])
        ax[1].set_title("Cell volume")
        ax[1].set_xlabel("Iteration")
        ax[1].set_ylabel("Å$^3$")
        ax[1].grid(True, alpha=0.3)
    else:
        ax[1].set_title("Cell volume (not recorded)")
        ax[1].axis("off")

    # ========== 3) |grad_frac| ==========
    if len(force_norm_vals) == len(iters):
        ax[2].plot(iters, force_norm_vals, "g-o", ms=3)
        add_anneal_lines(ax[2])
        ax[2].set_title("|grad_frac| (per-atom)")
        ax[2].set_xlabel("Iteration")
        ax[2].set_ylabel("Norm")
        ax[2].set_yscale("log")
        ax[2].grid(True, which="both", alpha=0.3)
    else:
        ax[2].set_title("|grad_frac| (not recorded)")
        ax[2].axis("off")

    # ========== 4) Energy Above Hull ==========
    ax[3].plot(iters, energy_above_hull_vals, "c-o", ms=3)
    ax[3].axhline(0.0, ls="--", color="gray", alpha=0.5)
    add_anneal_lines(ax[3])
    ax[3].set_title("Energy above hull (eV/atom)")
    ax[3].set_xlabel("Iteration")
    ax[3].set_ylabel("eV/atom")
    ax[3].grid(True, alpha=0.3)

    # ========== 5) Entropy ==========
    ax[4].plot(iters, ent_vals, "m-o", ms=3)
    add_anneal_lines(ax[4])
    ax[4].set_title("Entropy per atom")
    ax[4].set_xlabel("Iteration")
    ax[4].set_ylabel("Entropy/atom")
    ax[4].grid(True, alpha=0.3)

    # ========== 6) Final structure ==========
    ax_struct = ax[5]
    try:
        eat_final_np = Xf.detach().cpu().numpy()
        element_names = selected_elements

        elem_indices_vis = np.argmax(eat_final_np, axis=1)
        symbols_vis = [element_names[i] for i in elem_indices_vis]

        frac_vis = np.asarray(frac_final_np, dtype=float)
        cell_vis = np.asarray(cell_final, dtype=float)

        atoms_vis = Atoms(
            symbols=symbols_vis,
            scaled_positions=frac_vis,
            cell=cell_vis,
            pbc=True
        )

        unique_elems = list(dict.fromkeys(symbols_vis))
        color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
        elem_color_map = {
            elem: color_cycle[idx % len(color_cycle)] if color_cycle else "gray"
            for idx, elem in enumerate(unique_elems)
        }
        atom_colors = [elem_color_map[sym] for sym in symbols_vis]

        if ASE_PLOT_OK:
            plot_atoms(
                atoms_vis,
                ax_struct,
                radii=0.32,
                rotation="30x,25y,0z",
                colors=atom_colors,
            )
            ax_struct.set_axis_off()
            ax_struct.set_title("Final structure (argmax of X)")

            patches = [mpatches.Patch(color=elem_color_map[e], label=e) for e in unique_elems]
            star_patch = mpatches.Patch(color="k", label="⭐ Relaxed state")
            ax_struct.legend(handles=patches + [star_patch], loc="upper right", fontsize=8, frameon=True)
        else:
            cart = frac_vis @ cell_vis
            colors = [elem_color_map[s] for s in symbols_vis]
            ax_struct.scatter(cart[:, 0], cart[:, 1], s=35, c=colors,
                              alpha=0.8, edgecolors="k", linewidths=0.4)
            ax_struct.set_axis_off()
            ax_struct.set_title("Final structure (scatter fallback)")

    except Exception as e:
        ax_struct.axis("off")
        ax_struct.set_title(f"Structure plot failed: {e}")

    # ========== 7,8) Spares ==========
    ax[6].axis("off")
    ax[7].axis("off")

    plt.tight_layout()
    plot_path = log_dir / "optimization_dashboard.png"
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    log.info(f"📈 Saved optimization dashboard with ⭐ relaxed overlay to {plot_path}")

else:
    log.info("No history recorded; skipping plot generation.")
