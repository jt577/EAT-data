import numpy as np
import torch
from modules.minimizer import minimizer_EAT
import modules.logging as logging
from modules.run_dft import run_dft

class syntropizer_EAT:
    """
    syntropizer class for syntropizing (optimizing on discrete feasible set) a given function using Tsallis entropy.
    Used for optimizing DFT directly.
    Methods:
        _entropy(x): Computes the Tsallis q-Entropy and its gradient.
        _objfun(x): Computes the objective function and its gradient.
        _project(x): Projects x onto simplexes.
        _callback(x, f, err, k, converged): Callback function for optimization.
        _objfun_min(): Minimizes the quadratic form and performs entropy annealing.
    """

    def __init__(self, func, N, S, n_alpha, log_file, out_dir, dft_names, dft_infiles, dft_parallel, q, eta_init, n_anneals, anneal_mult, BFGS_hist, n_iter, tol, max_line_search, jitter, x0=None):
        """
        Attributes:
            func (callable): The DFT energy function to maximize.
            N (int): Number of atoms.
            S (int): Number of species.
            n_alpha (list): List of number of atoms we constrain per species.
            log_file (str): Path to log file.
            out_dir (str): Directory to save DFT outputs.
            dft_names (list): List of DFT code names.
            dft_infiles (list): List of DFT input file paths.
            dft_parallel (int): Number of parallel threads in DFT (for if you have > 1 k-point).
            n_vertices (int): Number of vertices to find.
            q (float): Tsallis entropy parameter.
            eta_init (float): Initial value for entropy annealing.
            eta (float): Current value for entropy annealing.
            n_anneals (int): Number of annealing steps.
            anneal_mult (float): Multiplicative factor to increase eta after each annealing step.
            BFGS_hist (int): History size for BFGS.
            n_iter (int): Maximum number of iterations for BFGS.
            tol (float): Tolerance for convergence in BFGS.
            max_line_search (int): Maximum number of line search steps.
            jitter (float): Jitter to add to initial x.
        """
        self.func = func
        self.N = N
        self.S = S
        self.n_alpha = n_alpha
        self.log_file = log_file
        self.out_dir = out_dir
        self.dft_names_original = dft_names
        self.dft_infiles = dft_infiles
        self.dft_parallel = dft_parallel
        self.q = q
        self.eta_init = eta_init
        self.eta = self.eta_init
        self.n_anneals = n_anneals
        self.anneal_mult = anneal_mult
        self.BFGS_hist = BFGS_hist
        self.n_iter = n_iter
        self.tol = tol
        self.max_line_search = max_line_search
        self.jitter = jitter
        self.x0 = x0

    def _run_DFT(self, x):
        """
        Objective for ordinary EAT optimization within DFT directly.
        """
        x_full = x.view(self.N, self.S)
        PBFGS_iter = self.function_minimizer.PBFGS_iteration
        try:
            line_search_iter = self.function_minimizer.line_search_iteration
        except:
            line_search_iter = 0
        run_name = f'EAT_anneal_{self.anneal_iter}_PBFGS_{PBFGS_iter}_lineSearch_{line_search_iter}'
        dft = run_dft(self.out_dir, run_name, self.dft_names_original, self.dft_infiles, self.dft_parallel, x_full)
        dft.run()
        self.dft_names = dft.dft_names  # Update DFT names in case they changed
        energies, energy_grads = dft.read()
        # Only return if DFT converged (sometimes it doesn't for some weird reason and this will mess up opt)
        if energies is not None and energy_grads is not None:
            y, grad_y = self.func(x_full, energies, energy_grads, self.dft_names)
            return y, grad_y.flatten()
        else:
            logging.log(f"     Warning: DFT did not converge for this x, returning None.", log_file=self.log_file)
            return None, None

    def _entropy(self, x):
        """
        Computes the Tsallis q-Entropy and its gradient.
        """
        x_full = x.view(self.N, -1)
        t = 0
        grad_t = torch.zeros_like(x_full)
        for i in range(self.N):
            if self.q==1:
                t += -torch.sum(x_full[i, :] * torch.log(x_full[i, :]+1e-10))
                grad_t[i, :] = -(torch.log(x_full[i, :] + 1e-10) + 1)
            else: 
                t += 1/(self.q-1) * (1-torch.sum(x_full[i, :]**self.q))
                grad_t[i, :] = -self.q/(self.q-1) * x_full[i, :]**(self.q-1)

        return t, grad_t.flatten()

    def _objfun(self, x):
        f, df = self._run_DFT(x)
        e, de = self._entropy(x)
        cost = -torch.tensor(f, dtype=torch.float64) + self.eta * e
        dcost = -torch.tensor(df , dtype=torch.float64) + self.eta * de
        return cost, dcost

    # ───────────────────────── helpers: simplex projections ─────────────────────────
    @staticmethod
    def _project_simplex_rows(M):
        """
        Project each row of M onto {x>=0, sum x = 1}. Returns a new array.
        """
        # M: (N, S)
        N, S = M.shape
        U = np.sort(M, axis=1)[:, ::-1]
        cssv = np.cumsum(U, axis=1) - 1
        ind = (np.arange(S) + 1)[None, :]  # shape (1,S)
        cond = U - cssv / ind > 0
        rho = cond.argmax(axis=1)  # first True along reversed? We'll compute last True robustly:
        # More robust: get last True index per row
        rho = cond.shape[1] - 1 - np.flip(cond, axis=1).argmax(axis=1)
        theta = cssv[np.arange(N), rho] / (rho + 1)
        return np.maximum(M - theta[:, None], 0.0)

    @staticmethod
    def _project_simplex_cols(M, sums):
        """
        Project each column j of M onto {x>=0, sum x = sums[j]}. Returns a new array.
        """
        N, S = M.shape
        out = np.empty_like(M)
        for j in range(S):
            v = M[:, j]
            s = float(sums[j])
            # standard simplex projection with target sum s
            u = np.sort(v)[::-1]
            cssv = np.cumsum(u) - s
            ind = np.arange(1, N+1)
            cond = u - cssv / ind > 0
            rho = cond.shape[0] - 1 - cond[::-1].argmax()  # last True
            theta = cssv[rho] / (rho + 1)
            out[:, j] = np.maximum(v - theta, 0.0)
        return out

    def _project(self, x):
        """
        Project onto:
          A: row-simplex constraints (sum over species = 1, nonneg)
          B: column-simplex constraints (sum over atoms = n_alpha, nonneg)
        If self.n_alpha is None, fall back to row-simplex only (old behavior).
        """
        x_np = x.detach().cpu().numpy().astype(np.float64, copy=False)
        M = x_np.reshape(self.N, self.S)

        if self.n_alpha is None:
            # old fast per-row simplex projection
            M_proj = self._project_simplex_rows(M)
            return torch.tensor(M_proj.reshape(-1), dtype=torch.float64)

        # Dykstra's algorithm for two convex sets A (rows) and B (cols)
        A_proj = self._project_simplex_rows
        def B_proj(Z): return self._project_simplex_cols(Z, self.n_alpha)

        Y = M.copy()
        P = np.zeros_like(M)  # correction for A
        Q = np.zeros_like(M)  # correction for B

        max_iters = 1000
        tol = 1e-10
        prev = Y.copy()

        for k in range(max_iters):
            # Project onto A with correction P
            XA = A_proj(Y + P)
            P = (Y + P) - XA

            # Project onto B with correction Q
            Y_new = B_proj(XA + Q)
            Q = (XA + Q) - Y_new

            # check convergence (relative Fro norm)
            diff = np.linalg.norm(Y_new - prev, ord='fro')
            base = max(1.0, np.linalg.norm(prev, ord='fro'))
            if diff / base < tol:
                Y = Y_new
                break

            prev = Y_new
            Y = Y_new

        # (Optional) final tiny clean-up to enforce sums numerically
        # Rows:
        row_sums = Y.sum(axis=1, keepdims=True)
        Y = np.maximum(Y, 0.0)
        Y = Y / np.clip(row_sums, 1e-16, None)  # rows to ~1
        # Columns: distribute tiny residual to hit n_alpha exactly
        col_sums = Y.sum(axis=0, keepdims=True)
        scale = np.clip(self.n_alpha / np.clip(col_sums, 1e-16, None), 0.0, np.inf)
        Y = Y * scale  # rescales each column; row sums may drift by ~1e-12

        return torch.tensor(Y.reshape(-1), dtype=torch.float64)

    def objfun_min(self):
        """
        Minimize and entropy anneal.
        """
        if self.x0 is not None:
            # Initialize at provided starting point with random jitter
            x_full = self.x0.view(self.N, self.S) + self.jitter * torch.randn(self.N, self.S, dtype=torch.float64)
        else:
            # Initialize at maximum entropy with random jitter
            x_full = 1/self.S * torch.ones(self.N, self.S) + self.jitter * torch.randn(self.N, self.S, dtype=torch.float64)
        x = x_full.flatten()
        # Project x
        x = self._project(x)

        eta = self.eta
        for anneal_iter in range(self.n_anneals):
            self.eta = eta
            self.anneal_iter = anneal_iter+1
            logging.log(f"\n\n=== Annealing step {self.anneal_iter}, eta = {self.eta:.6f} ===", log_file=self.log_file)
            self.function_minimizer = minimizer_EAT(self._objfun, x, args=(), proj_args=(), project=self._project, 
                                            log_file=self.log_file, n_atoms=self.N, n_species=self.S, q=self.q, eta=self.eta,
                                            BFGS_hist=self.BFGS_hist, maxit=self.n_iter, tol=self.tol, max_line_search=self.max_line_search)
            f_val, x = self.function_minimizer.PBFGS()
            e = self._entropy(x)[0]
            if e < 1e-3:
                logging.log("     Entropy is sufficiently low, optimization finished.", log_file=self.log_file)
                return f_val, x
                break
            logging.log(f"\nEnd of annealing step {self.anneal_iter}.\n", log_file=self.log_file)
            eta *= self.anneal_mult
        logging.log("\nFinal results", log_file=self.log_file)
        entropy, _ = self._entropy(x)
        objective = -f_val + self.eta * entropy
        logging.log(f"     f (minimizing): {f_val:.4e}", log_file=self.log_file)
        logging.log(f"     Objective (maximizing) = {objective:.4e}, Entropy = {entropy:.4e}.", log_file=self.log_file)
        logging.log("     x:", log_file=self.log_file)
        x_full = x.reshape(self.N, self.S)
        for row in x_full:
            logging.log("        " + "  ".join([f"{elem:.4f}" for elem in row]), log_file=self.log_file)
        return f_val, x