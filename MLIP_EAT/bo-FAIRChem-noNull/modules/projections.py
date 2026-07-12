import numpy as np
import torch

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


# modules/projections.py (or wherever project_flat_all lives)
def project_flat_all(z, n_x_atoms, n_elements, uc_lower, uc_upper, fix):
    n_x = n_x_atoms * n_elements
    if fix:
        x = z[:n_x].view(n_x_atoms, n_elements)
        x = project_simplex_rows(x).view(-1)
        return x
    else:
        x = z[:n_x].view(n_x_atoms, n_elements)
        x = project_simplex_rows(x).view(-1)
        uc = z[n_x:n_x+6]
        uc = torch.max(torch.min(uc, uc_upper), uc_lower)
        f = z[n_x+6:]
        return torch.cat([x, uc, f], dim=0)
