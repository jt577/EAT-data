import torch

# ============================================================================
# SIMPLEX PROJECTION
# ============================================================================

class projector:
    def __init__(self, n_x_atoms, n_elements, fix, mu=None, uc_lower=None, uc_upper=None):
        self.n_x_atoms = n_x_atoms
        self.n_elements = n_elements
        self.fix = fix
        self.mu = mu
        self.uc_lower = uc_lower
        self.uc_upper = uc_upper
    
    def project(self, X):
        """
        Apply alternating projections to the input tensor X, if mu is not None.
        """
        n_x = self.n_x_atoms * self.n_elements
        if self.fix:
            x = X[:n_x].view(self.n_x_atoms, self.n_elements)
            x = self._project_simplex_rows(x).view(-1)
            return x
        else:
            x = X[:n_x].view(self.n_x_atoms, self.n_elements)
            x = self._project_simplex_rows(x).view(-1)
            uc = X[n_x:n_x+6]
            uc = torch.max(torch.min(uc, self.uc_upper), self.uc_lower)
            f = X[n_x+6:]
            return torch.cat([x, uc, f], dim=0)

    def _project_simplex_rows(self, X):
        """
        Row-wise Euclidean projection onto simplex. Shape of input should be n x m, where m > 1.
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


















# DYKSTRA PROJECTION FOR ROWS + COLS SIMPLEX
# import torch

# class projector:
#     def __init__(self, n_x_atoms, n_elements, fix, mu=None, uc_lower=None, uc_upper=None):
#         self.n_x_atoms = n_x_atoms
#         self.n_elements = n_elements
#         self.fix = fix
#         self.mu = mu
#         self.uc_lower = uc_lower
#         self.uc_upper = uc_upper

#     def project(self, X: torch.Tensor) -> torch.Tensor:
#         n_x = self.n_x_atoms * self.n_elements

#         x = X[:n_x].view(self.n_x_atoms, self.n_elements)

#         if self.mu is not None:
#             x = self._dykstra_rows_cols(x, self.mu)

#         x = x.reshape(-1)

#         if self.fix:
#             return x

#         uc = X[n_x:n_x+6]
#         uc = torch.max(torch.min(uc, self.uc_upper), self.uc_lower)
#         f = X[n_x+6:]
#         return torch.cat([x, uc, f], dim=0)

#     # -----------------------------
#     # Dykstra for rows+cols simplex
#     # -----------------------------
#     def _dykstra_rows_cols(
#         self,
#         M: torch.Tensor,      # (N, S) = (n_x_atoms, n_elements)
#         mu: torch.Tensor,     # (S,)
#         max_iters: int = 2000,
#         tol: float | None = None,
#     ) -> torch.Tensor:
#         """
#         Project onto intersection:
#           A: rows sum to 1, nonneg
#           B: cols sum to n_x_atoms*mu, nonneg
#         using Dykstra's algorithm (Euclidean projections).
#         """
#         dtype = M.dtype
#         device = M.device

#         mu = mu.to(device=device, dtype=dtype)
#         mu = torch.clamp(mu, min=0.0)
#         mu_sum = mu.sum()
#         if mu_sum <= 0:
#             raise ValueError("mu must have positive sum.")
#         mu = mu / mu_sum  # IMPORTANT: make targets consistent

#         col_targets = self.n_x_atoms * mu  # (S,)
#         # sum(col_targets) == n_x_atoms required for feasibility

#         if tol is None:
#             tol = 1e-7 if dtype == torch.float32 else 1e-10

#         Y = M
#         P = torch.zeros_like(M)  # correction for A
#         Q = torch.zeros_like(M)  # correction for B

#         prev = Y
#         for _ in range(max_iters):
#             # A-projection (rows -> sum 1)
#             XA = self._project_simplex_rows_sum(Y + P, 1.0)
#             P = (Y + P) - XA

#             # B-projection (cols -> col_targets)
#             Y_new = self._project_simplex_cols_sum(XA + Q, col_targets)
#             Q = (XA + Q) - Y_new

#             # relative Frobenius change
#             diff = torch.linalg.norm(Y_new - prev)
#             base = torch.maximum(torch.tensor(1.0, device=device, dtype=dtype), torch.linalg.norm(prev))
#             if (diff / base) < tol:
#                 Y = Y_new
#                 break

#             prev = Y_new
#             Y = Y_new

#         return Y

#     # -----------------------------
#     # Simplex projections
#     # -----------------------------
#     @staticmethod
#     def _project_simplex_rows_sum(M: torch.Tensor, z) -> torch.Tensor:
#         """
#         Project each row i of M onto {x>=0, sum x = z_i}.
#         z can be scalar or shape (N,).
#         """
#         N, S = M.shape
#         device, dtype = M.device, M.dtype

#         if not torch.is_tensor(z):
#             z = torch.full((N,), float(z), device=device, dtype=dtype)
#         else:
#             z = z.to(device=device, dtype=dtype).view(-1)
#             if z.numel() != N:
#                 raise ValueError(f"z must be scalar or have shape ({N},)")

#         # Sort rows descending
#         U, _ = torch.sort(M, dim=1, descending=True)
#         cssv = torch.cumsum(U, dim=1) - z.view(-1, 1)

#         ind = torch.arange(1, S + 1, device=device, dtype=dtype).view(1, -1)
#         cond = U - cssv / ind > 0

#         # k = last index where cond is True, per row
#         # (flip then argmax gives first True in reversed order)
#         cond_flip = torch.flip(cond, dims=[1])
#         k_last = (S - 1) - torch.argmax(cond_flip.to(torch.int64), dim=1)
#         # If a row has no True (shouldn't happen), clamp
#         k_last = torch.clamp(k_last, min=0)

#         theta = cssv[torch.arange(N, device=device), k_last] / (k_last.to(dtype) + 1.0)
#         return torch.clamp(M - theta.view(-1, 1), min=0.0)

#     @staticmethod
#     def _project_simplex_cols_sum(M: torch.Tensor, sums: torch.Tensor) -> torch.Tensor:
#         """
#         Project each column j of M onto {x>=0, sum x = sums[j]}.
#         Uses row-simplex code on transpose (vectorized).
#         """
#         sums = sums.to(device=M.device, dtype=M.dtype).view(-1)  # (S,)
#         # project rows of M^T onto simplex with per-row sums = sums
#         MT = M.t()  # (S, N)
#         MT_proj = projector._project_simplex_rows_sum(MT, sums)  # each "row" is a column of M
#         return MT_proj.t()
