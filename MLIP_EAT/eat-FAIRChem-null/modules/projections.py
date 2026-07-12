import torch


class projector:
    def __init__(self, n_x_atoms, n_elements, fix, mu=None, uc_lower=None, uc_upper=None, mix_is_null=None):
        self.n_x_atoms = n_x_atoms
        self.n_elements = n_elements
        self.fix = fix
        self.mu = mu
        self.uc_lower = uc_lower
        self.uc_upper = uc_upper
        self.mix_is_null = mix_is_null

    def project(self, X):
        """
        Project variable rows:
          - normal EAT rows onto the simplex (sum = 1)
          - null rows onto the capped simplex (sum <= 1), with
            p_empty represented implicitly as 1 - sum(species probabilities)
        """
        n_x = self.n_x_atoms * self.n_elements
        x = X[:n_x].view(self.n_x_atoms, self.n_elements)
        x = self._project_rows(x).view(-1)

        if self.fix:
            return x

        uc = X[n_x:n_x + 6]
        uc = torch.max(torch.min(uc, self.uc_upper), self.uc_lower)
        f = X[n_x + 6:]
        return torch.cat([x, uc, f], dim=0)

    def _project_rows(self, X):
        if self.mix_is_null is None:
            return self._project_simplex_rows(X)

        null_mask = torch.as_tensor(self.mix_is_null, device=X.device, dtype=torch.bool).view(-1)
        if null_mask.numel() != X.shape[0]:
            raise ValueError("mix_is_null must have one entry per variable row")

        X_proj = X.clone()
        if torch.any(~null_mask):
            X_proj[~null_mask] = self._project_simplex_rows(X[~null_mask])
        if torch.any(null_mask):
            X_proj[null_mask] = self._project_capped_simplex_rows(X[null_mask])
        return X_proj

    def _project_simplex_rows(self, X):
        X_sorted, _ = torch.sort(X, dim=1, descending=True)
        cumsum = torch.cumsum(X_sorted, dim=1)
        rho = torch.arange(1, X.shape[1] + 1, device=X.device).view(1, -1)

        cond = X_sorted - (cumsum - 1) / rho > 0
        k = torch.clamp(cond.sum(dim=1), min=1)
        idx = (k - 1).long()

        theta = (cumsum[torch.arange(X.shape[0], device=X.device), idx] - 1) / k
        X_proj = torch.clamp(X - theta.view(-1, 1), min=0.0)
        return X_proj

    def _project_capped_simplex_rows(self, X):
        X_proj = torch.clamp(X, min=0.0)
        row_sums = X_proj.sum(dim=1)
        overfull = row_sums > 1.0
        if torch.any(overfull):
            X_proj[overfull] = self._project_simplex_rows(X[overfull])
        return X_proj
