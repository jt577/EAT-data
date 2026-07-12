#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import torch
import yaml  # pip install pyyaml

# FairChem
from fairchem.core import FAIRChemCalculator, pretrained_mlip

# BoTorch / GPyTorch
from botorch.acquisition.analytic import LogExpectedImprovement
from botorch.fit import fit_gpytorch_mll
from botorch.models.gp_regression_mixed import MixedSingleTaskGP
from botorch.optim.optimize import optimize_acqf_discrete
from botorch.utils.transforms import standardize
from gpytorch.mlls.exact_marginal_log_likelihood import ExactMarginalLogLikelihood

# Your modules
from modules.atoms_helpers import build_base_atoms
from modules.energy_above_hull import HullEnergyCalculator, PhaseRecord
from modules.logging_utils import setup_logging
from modules.plotting import plot_phase_diagram_from_db
from modules.seeding import save_rng_state, load_rng_state


# -----------------------------
# YAML helpers / IO
# -----------------------------
def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping/object. Got: {type(data)}")
    return data


def _subst_a(obj: Any, a: float) -> Any:
    """Recursively replace the string '{a}' with the float a inside lists/dicts."""
    if isinstance(obj, str):
        if obj.strip() == "{a}":
            return float(a)
        return obj
    if isinstance(obj, list):
        return [_subst_a(x, a) for x in obj]
    if isinstance(obj, dict):
        return {k: _subst_a(v, a) for k, v in obj.items()}
    return obj


# -----------------------------
# Chemistry helpers
# -----------------------------
def _composition_from_symbols(symbols: list[str]) -> dict[str, float]:
    from collections import Counter

    c = Counter(symbols)
    n = float(sum(c.values()))
    return {el: cnt / n for el, cnt in sorted(c.items())}


def _formation_energy_per_atom(
    energy_total_ev: float,
    symbols: list[str],
    reference_energies: dict[str, float],
) -> float:
    n = len(symbols)
    if n == 0:
        raise ValueError("Empty structure (n_atoms=0).")
    e_pa = energy_total_ev / float(n)
    comp = _composition_from_symbols(symbols)  # fractions
    ref_term = 0.0
    for el, frac in comp.items():
        if el not in reference_energies:
            raise ValueError(
                f"reference_energies.json missing element '{el}'. "
                f"Available keys (sample): {sorted(reference_energies.keys())[:30]}"
            )
        ref_term += frac * float(reference_energies[el])
    return e_pa - ref_term


def _counts_from_fractions(mu: Sequence[float], n_sites: int) -> list[int]:
    """
    Convert target fractions into integer counts that sum exactly to n_sites.
    Largest remainder (Hamilton) apportionment.
    """
    mu = np.asarray(mu, dtype=float)
    if mu.ndim != 1:
        raise ValueError(f"mu must be 1D, got shape {mu.shape}")
    if n_sites <= 0:
        raise ValueError("n_sites must be > 0")
    s = float(mu.sum())
    if s <= 0:
        raise ValueError("mu must sum to > 0")
    mu = mu / s

    raw = mu * n_sites
    base = np.floor(raw).astype(int)
    rem = raw - base

    missing = int(n_sites - base.sum())
    if missing > 0:
        idx = np.argsort(-rem)[:missing]
        base[idx] += 1
    elif missing < 0:
        idx = np.argsort(rem)[:(-missing)]
        base[idx] -= 1

    if base.sum() != n_sites or np.any(base < 0):
        raise ValueError(f"Bad apportionment: counts={base.tolist()} sum={base.sum()} n_sites={n_sites}")
    return base.tolist()


def _hash_symbols(symbols: list[str]) -> str:
    h = hashlib.sha1((" ".join(symbols)).encode("utf-8")).hexdigest()
    return h[:10]


def _label_from_fractions(elements: list[str], fracs: list[float]) -> str:
    return " ".join([f"{elements[i]}{fracs[i]:.3f}" for i in range(len(elements))])


# -----------------------------
# Plot helpers
# -----------------------------
def plot_formation_energy_histogram(
    eforms: list[float],
    *,
    out_path: Path,
    bins: Any = 40,
    hist_range: Any = None,
    title: str | None = None,
) -> Path:
    import matplotlib.pyplot as plt

    if len(eforms) == 0:
        return out_path

    data = np.array(eforms, dtype=float)
    rng = None
    if hist_range is not None:
        rng = (float(hist_range[0]), float(hist_range[1]))

    plt.figure()
    plt.hist(data, bins=bins, range=rng)
    plt.xlabel("Formation energy (eV/atom)")
    plt.ylabel("Count")
    if title is not None:
        plt.title(title)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    return out_path


def plot_best_trace(y: list[float], *, out_path: Path, title: str) -> Path:
    import matplotlib.pyplot as plt

    if len(y) == 0:
        return out_path
    yy = np.array(y, dtype=float)
    best = np.minimum.accumulate(yy)

    plt.figure()
    plt.plot(np.arange(1, len(best) + 1), best, linewidth=2.0)
    plt.xlabel("Evaluation #")
    plt.ylabel("Best formation energy so far (eV/atom)")
    plt.title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    return out_path


# -----------------------------
# Endmember injection
# -----------------------------
def ensure_endmembers_in_db_from_yaml(
    db_path: Path,
    design_space: list[str],
    endmembers: list[dict[str, Any]],
    *,
    overwrite: bool = False,
) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if not db_path.exists():
        db_path.write_text("{}", encoding="utf-8")

    db = json.loads(db_path.read_text(encoding="utf-8") or "{}")

    for em in endmembers:
        fractions = [float(x) for x in em["fractions"]]
        e_form = float(em["formation_energy_per_atom"])
        label = " ".join([f"{design_space[j]}{fractions[j]:.3f}" for j in range(len(design_space))])

        if (label in db) and (not overwrite):
            continue

        db[label] = {
            "elements": design_space,
            "fractions": fractions,
            "formation_energy_per_atom": e_form,
        }

    db_path.write_text(json.dumps(db, indent=2), encoding="utf-8")


# -----------------------------
# BO config
# -----------------------------
@dataclass(frozen=True)
class BOInputs:
    n_init: int
    n_total: int
    n_candidates: int
    n_neighbors: int
    neighbor_swaps: int
    seed: int



@dataclass(frozen=True)
class RunInputs:
    out_dir: Path
    compositions: list[list[float]]
    fix: bool
    selected_elements: list[str]
    lattice_vectors: list[list[float]]
    frac_coords_init: Any
    supercell_size: tuple[int, int, int]
    db_path: Path
    reference_json_path: Path
    endmembers: list[dict[str, Any]]

    model_name: str
    task_name: str

    # plotting knobs
    plot: bool
    plot_hist: bool
    plot_phase: bool
    plot_best: bool
    hist_bins: Any
    hist_range: Any
    update_every: int

    bo: BOInputs


def load_inputs(infile: Path) -> RunInputs:
    cfg = _load_yaml(infile)
    main_dir = Path(__file__).resolve().parent

    outdir = cfg.get("outdir", "out_botorch")
    out_dir = Path(outdir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    compositions = cfg["compositions"]
    if not isinstance(compositions, list) or not compositions:
        raise ValueError("config: compositions must be a non-empty list of lists")

    fix = bool(cfg.get("fix", True))
    selected_elements = list(cfg["selected_elements"])

    supercell_list = cfg["supercell_size"]
    if len(supercell_list) != 3:
        raise ValueError("config: supercell_size must have length 3")
    supercell_size = (int(supercell_list[0]), int(supercell_list[1]), int(supercell_list[2]))

    frac_coords_init = cfg["frac_coords_init"]
    lattice_vectors_raw = cfg["lattice_vectors"]

    if "lattice_a" in cfg:
        a = float(cfg["lattice_a"])
        lattice_vectors_raw = _subst_a(lattice_vectors_raw, a)

    lattice_vectors = [[float(x) for x in row] for row in lattice_vectors_raw]

    hull_db_filename = cfg.get("hull_db_filename", "hull_db.json")
    db_path = (out_dir / hull_db_filename).resolve()

    # reference energies
    ref_override = cfg.get("reference_json_path", None)
    if ref_override is None:
        reference_json_path = (main_dir / "data" / "reference_energies.json").resolve()
    else:
        p = Path(ref_override).expanduser()
        reference_json_path = (p if p.is_absolute() else (infile.parent / p)).resolve()

    if not reference_json_path.exists():
        raise FileNotFoundError(f"Reference energies JSON not found: {reference_json_path}")

    endmembers = cfg.get("endmembers", None)
    if endmembers is None or not isinstance(endmembers, list) or len(endmembers) == 0:
        raise ValueError("config: endmembers must be a non-empty list.")

    fairchem_cfg = cfg.get("fairchem", {}) or {}
    model_name = str(fairchem_cfg.get("model_name", "uma-s-1"))
    task_name = str(fairchem_cfg.get("task_name", "omat"))

    # plotting knobs
    plot = bool(cfg.get("plot", True))
    plot_phase = bool(cfg.get("plot_phase", True))
    plot_hist = bool(cfg.get("plot_hist", True))
    plot_best = bool(cfg.get("plot_best", True))
    hist_bins = cfg.get("hist_bins", 40)
    hist_range = cfg.get("hist_range", None)

    update_every = int(cfg.get("update_every", 20))
    if update_every <= 0:
        raise ValueError("update_every must be >= 1")

    # BO knobs
    bo_cfg = cfg.get("bo", {}) or {}
    seed = int(bo_cfg.get("seed", cfg.get("random_seed", 0)))

    n_init = int(bo_cfg.get("n_init", 50))
    n_total = int(bo_cfg.get("n_total", 300))
    if n_init <= 0 or n_total <= 0 or n_init > n_total:
        raise ValueError("bo.n_init and bo.n_total must be >0 and n_init <= n_total")

    q = int(bo_cfg.get("q", 1))
    if q <= 0:
        raise ValueError("bo.q must be >= 1")

    n_candidates = int(bo_cfg.get("n_candidates", 4096))
    if n_candidates <= 0:
        raise ValueError("bo.n_candidates must be > 0")

    n_neighbors = int(bo_cfg.get("n_neighbors", 512))
    neighbor_swaps = int(bo_cfg.get("neighbor_swaps", 4))
    if n_neighbors < 0 or neighbor_swaps < 1:
        raise ValueError("bo.n_neighbors must be >=0 and bo.neighbor_swaps >=1")

    mc_samples = int(bo_cfg.get("mc_samples", 256))
    if mc_samples <= 0:
        raise ValueError("bo.mc_samples must be > 0")

    return RunInputs(
        out_dir=out_dir,
        compositions=[[float(v) for v in comp] for comp in compositions],
        fix=fix,
        selected_elements=selected_elements,
        lattice_vectors=lattice_vectors,
        frac_coords_init=frac_coords_init,
        supercell_size=supercell_size,
        db_path=db_path,
        reference_json_path=reference_json_path,
        endmembers=endmembers,
        model_name=model_name,
        task_name=task_name,
        plot=plot,
        plot_hist=plot_hist,
        plot_phase=plot_phase,
        plot_best=plot_best,
        hist_bins=hist_bins,
        hist_range=hist_range,
        update_every=update_every,
        bo=BOInputs(
            n_init=n_init,
            n_total=n_total,
            n_candidates=n_candidates,
            n_neighbors=n_neighbors,
            neighbor_swaps=neighbor_swaps,
            seed=seed,
        ),
    )


# -----------------------------
# Repair operator on integer codes (enforce exact counts)
# -----------------------------
def repair_codes_to_counts(
    proposed: np.ndarray,
    *,
    k: int,
    target_counts: list[int],
    seed: int,
) -> np.ndarray:
    """
    proposed: (n_sites,) int in [0, k-1]
    target_counts: length k, sum == n_sites

    Deterministic repair by flipping random sites from overrepresented -> underrepresented.
    """
    proposed = np.asarray(proposed, dtype=int).copy()
    n = int(proposed.shape[0])
    if n == 0:
        return proposed
    if len(target_counts) != k:
        raise ValueError("target_counts length != k")
    if sum(target_counts) != n:
        raise ValueError("target_counts must sum to n_sites")

    h = hashlib.sha1((proposed.tobytes() + str(seed).encode("utf-8"))).hexdigest()
    rng = np.random.default_rng(int(h[:16], 16))

    counts = np.bincount(proposed, minlength=k).tolist()
    surplus = [counts[i] - target_counts[i] for i in range(k)]

    where = [np.where(proposed == i)[0].tolist() for i in range(k)]

    while True:
        donors = [i for i in range(k) if surplus[i] > 0]
        receivers = [i for i in range(k) if surplus[i] < 0]
        if not donors and not receivers:
            break
        if not donors or not receivers:
            break

        d = int(rng.choice(donors))
        r = int(rng.choice(receivers))

        if not where[d]:
            where = [np.where(proposed == i)[0].tolist() for i in range(k)]

        pos = int(rng.choice(where[d]))
        proposed[pos] = r

        where[d].remove(pos)
        where[r].append(pos)
        surplus[d] -= 1
        surplus[r] += 1

    final_counts = np.bincount(proposed, minlength=k).tolist()
    if final_counts != list(target_counts):
        raise RuntimeError(f"repair failed: got counts {final_counts} expected {target_counts}")

    return proposed


def codes_to_symbols(codes: np.ndarray, elements: list[str]) -> list[str]:
    return [elements[int(c)] for c in codes.tolist()]


def hash_codes(codes: np.ndarray) -> str:
    return hashlib.sha1(codes.tobytes()).hexdigest()[:10]


# -----------------------------
# Candidate generation
# -----------------------------
def sample_random_candidates(
    *,
    n: int,
    n_sites: int,
    k: int,
    target_counts: list[int],
    seed: int,
) -> np.ndarray:
    """
    Returns (n, n_sites) int codes in [0,k-1], each repaired to target_counts.
    """
    rng = np.random.default_rng(seed)
    out = np.empty((n, n_sites), dtype=int)
    for i in range(n):
        prop = rng.integers(low=0, high=k, size=(n_sites,), dtype=int)
        out[i] = repair_codes_to_counts(prop, k=k, target_counts=target_counts, seed=seed + i)
    return out


def neighbor_swap_candidates(
    *,
    x_best: np.ndarray,
    n: int,
    k: int,
    target_counts: list[int],
    n_swaps: int,
    seed: int,
) -> np.ndarray:
    """
    Build candidates by doing random swaps in x_best (composition preserved automatically),
    then (optional) repair for safety.
    """
    rng = np.random.default_rng(seed)
    x_best = np.asarray(x_best, dtype=int)
    n_sites = int(x_best.shape[0])

    out = np.empty((n, n_sites), dtype=int)
    for i in range(n):
        x = x_best.copy()
        for _ in range(n_swaps):
            # pick two positions with different species
            for _try in range(50):
                a = int(rng.integers(0, n_sites))
                b = int(rng.integers(0, n_sites))
                if x[a] != x[b]:
                    x[a], x[b] = x[b], x[a]
                    break
        x = repair_codes_to_counts(x, k=k, target_counts=target_counts, seed=seed + 10000 + i)
        out[i] = x
    return out


def unique_rows_int(X: np.ndarray) -> np.ndarray:
    # stable uniqueness via byte-views
    if X.size == 0:
        return X
    b = X.view(np.dtype((np.void, X.dtype.itemsize * X.shape[1])))
    _, idx = np.unique(b, return_index=True)
    return X[np.sort(idx)]


# -----------------------------
# Evaluator
# -----------------------------
class BOCore:
    def __init__(
        self,
        *,
        log,
        atoms_base,
        calc,
        reference_energies: dict[str, float],
        selected_elements: list[str],
        target_counts: list[int],
        hullcalc: HullEnergyCalculator,
        db_json_path: str,
        out_dir: Path,
        comp_idx: int,
        supercell_size: tuple[int, int, int],
        plot_phase: bool,
        plot_hist: bool,
        plot_best: bool,
        hist_bins: Any,
        hist_range: Any,
        update_every: int,
        device: torch.device,
        seed: int,
    ):
        self.log = log
        self.atoms_base = atoms_base
        self.calc = calc
        self.reference_energies = reference_energies
        self.selected_elements = selected_elements
        self.k = len(selected_elements)
        self.target_counts = target_counts
        self.hullcalc = hullcalc
        self.db_json_path = db_json_path
        self.out_dir = out_dir
        self.comp_idx = comp_idx
        self.supercell_size = supercell_size

        self.plot_phase = plot_phase
        self.plot_hist = plot_hist
        self.plot_best = plot_best
        self.hist_bins = hist_bins
        self.hist_range = hist_range
        self.update_every = update_every

        self.device = device
        self.seed = int(seed)

        self.X_list: list[np.ndarray] = []  # each (n_sites,)
        self.Y_list: list[float] = []
        self.cache: dict[str, float] = {}  # hash10 -> y

    def _update_plots(self) -> None:
        n_done = len(self.Y_list)
        if self.plot_hist:
            plot_formation_energy_histogram(
                self.Y_list,
                out_path=self.out_dir / f"formation_energy_hist_comp{self.comp_idx:02d}.png",
                bins=self.hist_bins,
                hist_range=self.hist_range,
                title=f"Formation energies (comp {self.comp_idx})  n={n_done}",
            )
            self.log.info("📊 Updated histogram.")

        if self.plot_best:
            plot_best_trace(
                self.Y_list,
                out_path=self.out_dir / f"best_trace_comp{self.comp_idx:02d}.png",
                title=f"Best-so-far (comp {self.comp_idx})  n={n_done}",
            )
            self.log.info("🏁 Updated best trace.")

        if self.plot_phase and len(self.selected_elements) <= 2:
            plot_phase_diagram_from_db(
                db_json_path=self.db_json_path,
                design_space=self.selected_elements,
                out_dir=self.out_dir,
                fname="phase_diagram.png",
                annotate=False,
            )
            self.log.info("📈 Updated phase diagram.")

    def evaluate_codes(self, codes: np.ndarray) -> float:
        codes = np.asarray(codes, dtype=int).copy()
        codes = repair_codes_to_counts(
            codes,
            k=self.k,
            target_counts=self.target_counts,
            seed=self.seed,
        )
        h = hash_codes(codes)
        if h in self.cache:
            return float(self.cache[h])

        symbols = codes_to_symbols(codes, self.selected_elements)

        atoms = self.atoms_base.copy()
        atoms.set_chemical_symbols(symbols)
        atoms.calc = self.calc

        energy_total_ev = float(atoms.get_potential_energy())
        e_form_pa = _formation_energy_per_atom(
            energy_total_ev=energy_total_ev,
            symbols=symbols,
            reference_energies=self.reference_energies,
        )

        # write to DB
        comp_sym = _composition_from_symbols(symbols)
        fr_list = [float(comp_sym.get(el, 0.0)) for el in self.selected_elements]
        base_label = _label_from_fractions(self.selected_elements, fr_list)
        label = f"{base_label} botorch-{h}"

        record = PhaseRecord(
            elements=self.selected_elements,
            fractions=fr_list,
            formation_energy_per_atom=float(e_form_pa),
            total_energy_ev=float(energy_total_ev),
            n_atoms=len(symbols),
            symbols=list(symbols),
            frac_coords=atoms.get_scaled_positions().tolist(),
            lattice_vectors=atoms.get_cell().array.tolist(),
            cell_pbc=[True, True, True],
            supercell_size=list(self.supercell_size),
            metadata={"source": "botorch", "hash": h},
        )
        self.hullcalc.add_phase_to_db(label=label, record=record)

        # cache + traces
        self.cache[h] = float(e_form_pa)
        self.X_list.append(codes)
        self.Y_list.append(float(e_form_pa))

        n_done = len(self.Y_list)
        if n_done == 1 or n_done % 10 == 0:
            self.log.info(f"🧠 Eval {n_done:5d}  Eform={e_form_pa:.6f} eV/atom  hash={h}")

        if n_done == 1 or (n_done % self.update_every == 0):
            self._update_plots()

        return float(e_form_pa)

    def best_so_far(self) -> tuple[float, np.ndarray]:
        y = np.array(self.Y_list, dtype=float)
        i = int(np.argmin(y))
        return float(y[i]), np.asarray(self.X_list[i], dtype=int)


# -----------------------------
# GP + qEI step
# -----------------------------
def fit_model_mixedgp(train_X: torch.Tensor, train_Y: torch.Tensor, cat_dims: list[int]) -> MixedSingleTaskGP:
    # BoTorch expects train_Y shape (n, 1)
    if train_Y.ndim == 1:
        train_Y = train_Y.unsqueeze(-1)
    model = MixedSingleTaskGP(
        train_X=train_X,
        train_Y=train_Y,
        cat_dims=cat_dims,
    )
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    return model


def select_next_via_logei(
    *,
    model: MixedSingleTaskGP,
    train_Y: torch.Tensor,
    choices: torch.Tensor,
    X_avoid: torch.Tensor | None,
) -> torch.Tensor:
    """
    Analytic LogEI is SINGLE-POINT (q=1) and does NOT use MC sampling.
    choices: (N, d) tensor
    returns: (1, d) tensor
    """
    # We are minimizing formation energy.
    # BoTorch acquisitions assume "maximize", so we model NEGATIVE energy.
    best_f = train_Y.max().item()  # best observed value in maximization space

    acqf = LogExpectedImprovement(model=model, best_f=best_f)

    X_best, _ = optimize_acqf_discrete(
        acq_function=acqf,
        q=1,
        choices=choices,
        unique=True,
        X_avoid=X_avoid,
    )
    return X_best



# -----------------------------
# Main
# -----------------------------
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, required=True, help="YAML config file")
    args = parser.parse_args()

    infile = Path(args.infile).expanduser().resolve()
    if not infile.exists():
        raise FileNotFoundError(f"infile not found: {infile}")

    inp = load_inputs(infile)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log = setup_logging(inp.out_dir)
    log.info(f"Output dir: {inp.out_dir}")
    log.info(f"Device: {device}")
    log.info(f"Infile: {infile}")
    log.info(f"Reference energies: {inp.reference_json_path}")
    log.info(
        f"BoTorch BO: n_init={inp.bo.n_init} n_total={inp.bo.n_total} "
        f"n_candidates={inp.bo.n_candidates} n_neighbors={inp.bo.n_neighbors}"
    )


    # reference energies
    reference_energies = json.loads(inp.reference_json_path.read_text(encoding="utf-8"))

    # keep RNG stable around FairChem init
    rng_state = save_rng_state()
    predictor = pretrained_mlip.get_predict_unit(inp.model_name)
    calc = FAIRChemCalculator(predictor, task_name=inp.task_name)
    load_rng_state(rng_state)

    # DB init
    inp.db_path.parent.mkdir(parents=True, exist_ok=True)
    if not inp.db_path.exists():
        inp.db_path.write_text("{}", encoding="utf-8")
    DB_JSON_PATH = str(inp.db_path)

    ensure_endmembers_in_db_from_yaml(inp.db_path, inp.selected_elements, inp.endmembers)

    # atoms template
    atoms_base = build_base_atoms(inp.frac_coords_init, inp.lattice_vectors, inp.supercell_size)
    n_sites = len(atoms_base)
    log.info(f"n_sites in supercell structure: {n_sites}")

    # global dtype for GP
    torch.set_default_dtype(torch.double)

    # loop compositions
    for comp_idx, composition in enumerate(inp.compositions):
        if len(composition) != len(inp.selected_elements):
            raise ValueError(
                f"Composition length {len(composition)} != number of elements {len(inp.selected_elements)} "
                f"for composition={composition} elements={inp.selected_elements}"
            )

        # target exact counts
        counts = _counts_from_fractions(composition, n_sites)
        log.info(
            f"[comp {comp_idx}] target={composition} integer_counts={dict(zip(inp.selected_elements, counts))}"
        )

        hullcalc = HullEnergyCalculator(DB_JSON_PATH, inp.selected_elements)

        bo = BOCore(
            log=log,
            atoms_base=atoms_base,
            calc=calc,
            reference_energies=reference_energies,
            selected_elements=inp.selected_elements,
            target_counts=counts,
            hullcalc=hullcalc,
            db_json_path=DB_JSON_PATH,
            out_dir=inp.out_dir,
            comp_idx=comp_idx,
            supercell_size=inp.supercell_size,
            plot_phase=inp.plot_phase,
            plot_hist=inp.plot_hist,
            plot_best=inp.plot_best,
            hist_bins=inp.hist_bins,
            hist_range=inp.hist_range,
            update_every=inp.update_every,
            device=device,
            seed=inp.bo.seed + 1000 * comp_idx,
        )

        k = len(inp.selected_elements)

        # --- Initialization: random repaired orderings ---
        init_codes = sample_random_candidates(
            n=inp.bo.n_init,
            n_sites=n_sites,
            k=k,
            target_counts=counts,
            seed=inp.bo.seed + 12345 + 1000 * comp_idx,
        )
        init_codes = unique_rows_int(init_codes)
        log.info(f"Init: evaluating {init_codes.shape[0]} unique repaired random configs.")
        for i in range(init_codes.shape[0]):
            bo.evaluate_codes(init_codes[i])

        # --- BO loop ---
        while len(bo.Y_list) < inp.bo.n_total:
            n_done = len(bo.Y_list)
            best_y, best_x = bo.best_so_far()
            log.info(f"🔎 BO step: n_done={n_done}/{inp.bo.n_total}  best={best_y:.6f} eV/atom")

            # candidate pool: random + neighbors around best
            n_rand = inp.bo.n_candidates
            X_rand = sample_random_candidates(
                n=n_rand,
                n_sites=n_sites,
                k=k,
                target_counts=counts,
                seed=inp.bo.seed + 999 + n_done + 1000 * comp_idx,
            )

            X_pool = X_rand
            if inp.bo.n_neighbors > 0:
                X_nb = neighbor_swap_candidates(
                    x_best=best_x,
                    n=inp.bo.n_neighbors,
                    k=k,
                    target_counts=counts,
                    n_swaps=inp.bo.neighbor_swaps,
                    seed=inp.bo.seed + 4242 + n_done + 1000 * comp_idx,
                )
                X_pool = np.vstack([X_pool, X_nb])
            X_pool = unique_rows_int(X_pool)

            # --- after you build X_pool and filter vs obs_hashes ---

            # Canonicalize observed points
            X_obs = np.stack([
                repair_codes_to_counts(x, k=k, target_counts=counts, seed=bo.seed)
                for x in bo.X_list
            ], axis=0)

            # Canonicalize pool the same way evaluation will
            X_pool = np.stack([
                repair_codes_to_counts(row, k=k, target_counts=counts, seed=bo.seed)
                for row in X_pool
            ], axis=0)
            X_pool = unique_rows_int(X_pool)

            # tensors for GP (IMPORTANT: use canonical X_obs)
            train_X = torch.tensor(X_obs, device=device, dtype=torch.double)
            train_Y = torch.tensor(np.array(bo.Y_list, dtype=float), device=device, dtype=torch.double)

            train_Y_model = -train_Y  # maximize
            cat_dims = list(range(n_sites))

            model = fit_model_mixedgp(train_X, train_Y_model, cat_dims=cat_dims).to(device)

            choices = torch.tensor(X_pool, device=device, dtype=torch.double)
            X_avoid = torch.tensor(X_obs, device=device, dtype=torch.double)

            X_next = select_next_via_logei(
                model=model,
                train_Y=train_Y_model,   # maximize space
                choices=choices,
                X_avoid=X_avoid,
            )

            # evaluate selected points
            X_next_np = X_next.detach().cpu().numpy().round().astype(int)  # q x d
            for j in range(X_next_np.shape[0]):
                bo.evaluate_codes(X_next_np[j])
                if len(bo.Y_list) >= inp.bo.n_total:
                    break

        # final plot pass
        try:
            bo._update_plots()
        except Exception as e:
            log.warning(f"Could not finalize plots: {e}")

        best_y, _best_x = bo.best_so_far()
        log.info(f"🏁 Done comp {comp_idx}. Best observed formation energy = {best_y:.6f} eV/atom")

    # optional final phase diagram
    if inp.plot and len(inp.selected_elements) <= 2:
        plot_path = plot_phase_diagram_from_db(
            db_json_path=DB_JSON_PATH,
            design_space=inp.selected_elements,
            out_dir=inp.out_dir,
            fname="phase_diagram.png",
            annotate=False,
        )
        log.info(f"📈 Saved phase diagram to {plot_path}")

    log.info("✅ Done.")


if __name__ == "__main__":
    main()
