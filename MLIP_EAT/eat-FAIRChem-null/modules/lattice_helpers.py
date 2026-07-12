import numpy as np
import torch
import math

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