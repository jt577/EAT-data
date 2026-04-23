import torch
import modules.logging as logging

class minimizer_EAT:
    """
    minimizer_EAT class implements optimization for minimizing functions for ordinary EAT using projected limited memory BFGS method.
    This class provides methods to perform optimization on a given function with its gradient, utilizing techniques such as 
    spectral projected gradient (SPG) and line minimization. The class is designed to handle constraints by projecting 
    the iterates back onto the feasible set.
    Methods:
        PBFGS(fun_jac, x0, BFGS_hist, args, proj_args, callback, project, maxit=100, tol=1e-10):
            Performs projected limited memory BFGS optimization.
        _construct_B(s, y, S_list, Y_list):
            Constructs the matrices needed to obtain the limited memory BFGS matrix.
        _construct_r(delta, N, M_inv, s, y):
            Constructs the damped BFGS update version of the gradient difference.
        _SPG(x0, f0, g0, delta, N, M_inv, SPG_hist, maxit=1000, tol=1e-3):
            Performs spectral projected gradient optimization.
        _line_min(fun_jac, x0, f0, g0, p, f_b, args, amax=1, c1=0.0001, maxiter=10):
            Performs line minimization to find the optimal step size.
        _cubic_interp(an, fn, gn, ao, fo, go):
            Obtains a new step size by cubic interpolation.
    """

    def __init__(self, fun_jac, x0, args, project, proj_args, log_file, n_atoms, n_species, q, eta, BFGS_hist=50, maxit=100, tol=1e-10, max_line_search=10):
        """
        Initializes the minimizer class.

        Attributes:
            fun_jac: Function to minimize (along with its gradient)
            x0: Initial guess
            args: Extra args passed to fun and jac
            project: Projection back onto nearest point on valid domain
            proj_args: Arguments for the projection function
            callback: Callback function for displaying information at each iteration
            log_file: File to log output for ordinary EAT
            n_atoms: Number of atoms
            n_species: Number of species
            q: Tsallis entropy parameter
            eta: Entropy annealing parameter
            BFGS_hist: History length for function values
            maxit: Maximum iterations
            tol: Tolerance of convergence
            max_line_search: Maximum line search iterations
        """
        self.fun_jac = fun_jac
        self.x0 = x0
        self.args = args
        self.proj_args = proj_args
        self.log_file = log_file
        self.n_atoms = n_atoms
        self.n_species = n_species
        self.q = q
        self.eta = eta
        self.project = project
        self.BFGS_hist = BFGS_hist
        self.maxit = maxit
        self.tol = tol
        self.max_line_search = max_line_search
    
    def PBFGS(self):
        """
        Perform projected limited memory BFGS.
        """
    
        device, dtype = self.x0.device, self.x0.dtype
        x = self.x0.clone()
        f_new = None
        g_new = None
        k = 1
        self.PBFGS_iteration = 0
        S_list = []
        R_list = []
        converged = False
        reset = True  # flag to indicate if we need to reset the quasi-Hessian (or begin optimization)
        while k <= self.maxit:
            if self.PBFGS_iteration != k: # only update if new iteration (avoids double printing PBFGS iteration headline in output)
                self.PBFGS_iteration = k
                # log iteration to log file if doing ordinary EAT
                logging.log(f"\nPBFGS iteration {k}", log_file=self.log_file)

            if reset == True:
                # Perform gradient descent step on first iteration to initialize delta
                f, g = self.fun_jac(x, *self.args)
                err = torch.norm(self.project(x - g, *self.proj_args) - x).item()
                if isinstance(f, torch.Tensor):
                    f = f.item()
                alpha = 1e-3  # initial step size for gradient descent (TODO: optimize this automatically instead of hardcoding)
                x = x - alpha * g
                x = self.project(x, *self.proj_args)  # project again to make sure no negative values
                f_new, g_new = self.fun_jac(x, *self.args)
                entropy = self._entropy(x)
                objective = -f_new + self.eta * entropy
                logging.log(f"     Initial gradient descent step.", log_file=self.log_file)
                logging.log(f"     f (minimizing) = {f_new:.4e}, |P(x-g)-x| = {err:.3e}.", log_file=self.log_file)
                logging.log(f"     Objective (maximizing) = {objective:.4e}, Entropy = {entropy:.4e}.", log_file=self.log_file)
                logging.log(f"     x:", log_file=self.log_file)
                for row in x.view(self.n_atoms, self.n_species):
                    logging.log("        " + "  ".join([f"{elem:.3f}" for elem in row]), log_file=self.log_file)
                s = -alpha * g
                y = g_new - g
                r = y.clone()  # use y as initial r
                # Set delta, N, M_inv for the first iteration
                delta = (y.dot(y) / y.dot(s)).item() if y.dot(s).item() > 0 else 1.0
                N = torch.zeros((len(x), 1), device=device, dtype=dtype)
                M_inv = torch.zeros((1, 1), device=device, dtype=dtype)
                # set reset to false
                reset = False
            else:
                f = f_new
                g = g_new

                err = torch.norm(self.project(x - g, *self.proj_args) - x).item()
                if err < self.tol:
                    # converged
                    converged = True
                    break
                
                x_star, SPG_converged = self._SPG(
                    x0=x, f0=f, g0=g, delta=delta,
                    N=N, M_inv=M_inv, SPG_hist=3 # hard-code mild SPG history length
                )

                if SPG_converged == False:
                    # SPG did not converge likely due to poor conditioning; reset quasi-Hessian
                    logging.log('     Warning: SPG did not converge, resetting quasi-Hessian.', log_file=self.log_file)
                    S_list = []
                    R_list = []
                    reset = True
                    continue
                else:
                    x_star = self.project(x_star, *self.proj_args)  # project result to be safe
                    d = x_star - x

                    alpha, f_temp, g_temp = self._line_min(
                        fun_jac=self.fun_jac, x0=x, f0=f, g0=g, p=d,
                        f_b=f, args=self.args, amax=1, c1=1e-4, 
                        maxiter=self.max_line_search, verbal=True,
                    )
                    if (alpha is None) or (alpha < 0) or isinstance(alpha, complex):        # Line search failed.
                        # If ordinary EAT then quit optimization and increase entropy annealing parameter (eta)
                        logging.log(f"     Warning: Line search failed, likely due to bad conditioning. Increasing entropy parameter eta.", log_file=self.log_file)
                        break

                    g_new = g_temp
                    f_new = f_temp
                    x = x + alpha * d
                    x = self.project(x, *self.proj_args)  # project again to make sure no negative values
                    entropy = self._entropy(x)
                    objective = -f_new + self.eta * entropy
                    logging.log(f"     f (minimizing) = {f_new:.4e}, |P(x-g)-x| = {err:.3e}.", log_file=self.log_file)
                    logging.log(f"     Objective (maximizing) = {objective:.4e}, Entropy = {entropy:.4e}.", log_file=self.log_file)
                    logging.log(f"     x:", log_file=self.log_file)
                    for row in x.view(self.n_atoms, self.n_species):
                        logging.log("        " + "  ".join([f"{elem:.3f}" for elem in row]), log_file=self.log_file)

                    s = alpha * d
                    y = g_new - g
                    r = self._construct_r(delta, N, M_inv, s, y) # construct r for damped update

                    # Construct BFGS matrix
                    delta, N, M_inv, S_list, R_list = self._construct_B(s, r, S_list, R_list)

            # After updating x, logging, and BFGS matrices:
            f = f_new
            g = g_new
            k += 1

        return f, x

    def _construct_B(self, s, y, S_list, Y_list):
        """
        Construct the matrices needed to obtain the limited memory BFGS matrix.

        Args:
            s: Current step (descent direction)
            y: Gradient difference
            S_list: List of previous steps
            Y_list: List of previous gradient differences
        Returns:
            delta, N, M_inv: Matrices used to update the BFGS matrix by B = delta - N @ inv(M) @ N.T
        """
        device = s.device
        dtype = s.dtype
        if len(S_list) < self.BFGS_hist:
            S_list.append(s.clone())
            Y_list.append(y.clone())
        else:
            S_list.pop(0)
            S_list.append(s.clone())
            Y_list.pop(0)
            Y_list.append(y.clone())

        S = torch.stack(S_list, dim=1)
        Y = torch.stack(Y_list, dim=1)
        gamma = torch.dot(s, y).item() / torch.dot(y, y).item()
        delta = 1.0 / gamma

        L = torch.zeros((S.shape[1], S.shape[1]), device=device, dtype=dtype)
        for i in range(S.shape[1]):
            for j in range(S.shape[1]):
                if i > j:
                    L[i, j] = torch.dot(S_list[i], Y_list[j]).item()

        D = torch.zeros((S.shape[1], S.shape[1]), device=device, dtype=dtype)
        for i in range(S.shape[1]):
            D[i, i] = torch.dot(S_list[i], Y_list[i]).item()

        N = torch.cat([delta * S, Y], dim=1)
        A = delta * S.t().matmul(S)
        top = torch.cat([A, L], dim=1)
        bottom = torch.cat([L.t(), -D], dim=1)
        M = torch.cat([top, bottom], dim=0)
        M_inv = torch.inverse(M)
        return delta, N, M_inv, S_list, Y_list

    def _construct_r(self, delta, N, M_inv, s, y):
        """
        Construct the damped BFGS update version of y.

        Args:
            delta, N, M_inv: Necessary to build quasi-Hessian
            s: Current step (descent direction)
            y: Gradient difference
        """
        # Compute B @ s
        t = N.t() @ s
        u = M_inv @ t            
        v = N @ u                
        Bs = delta * s - v       

        # Compute theta based on the curvature condition
        if s @ y >= 0.2 * s @ Bs:
            theta = 1
        else:
            theta = 0.8 * s @ Bs / (s @ Bs - s @ y)

        # Compute r
        r = theta * y + (1 - theta) * Bs

        return r

    def _SPG(self, x0, f0, g0, delta, N, M_inv, SPG_hist, maxit=1000, tol=1e-3):
        """
        Spectral Projected Gradient for ½‖x-x0‖_H² + g0ᵀ(x-x0) + f0.

        Args:
            x0: Initial point
            f0: Initial function value at x0
            g0: Gradient at x0
            delta, N, M_inv: Matrices for limited memory BFGS update
            c: Inverse volume fraction of the design domain
            SPG_hist: History length for function values
            maxit: Maximum number of iterations
            tol: Tolerance for convergence
        """
        def quad(x):
            # Compute the quadratic function value and gradient
            r = N.t().matmul(x - x0)
            r = M_inv.matmul(r)
            r = N.matmul(r)
            r = delta * (x - x0) - r
            return (0.5 * torch.dot(x - x0, r) + torch.dot(g0, x - x0)).item() + f0, r + g0

        a_max = 1e10
        a_min = 1e-10
        x = self.project(x0, *self.proj_args)           # ensure feasibility
        f, g = quad(x)              # initial function value
        a_bb = 1.0
        f_hist = [f]  # history of function values

        for k in range(maxit):
            a_k = min(a_max, max(a_min, a_bb))
            d_k = self.project(x - a_k * g, *self.proj_args) - x  # descent direction
            f_b = max(f_hist)
            a, f_new_val, g_new = self._line_min(
                fun_jac=quad, x0=x, f0=f, g0=g,
                p=d_k, f_b=f_b, args=(), amax=1, c1=1e-4, 
                maxiter=10, verbal=False,
            )
            # if line search failed, use a = 1 (fall back to bb step since it is encoded already in d_k)
            if (a is None) or (a < 0) or isinstance(a, complex):
                a = 1.0
                f_new_val, g_new = quad(x + a * d_k)  # re-evaluate function and gradient

            if len(f_hist) < SPG_hist:
                f_hist.append(f_new_val)
            else:
                f_hist = f_hist[1:] + [f_new_val]

            x_new = x + a * d_k
            s = x_new - x
            y = g_new - g
            num = torch.dot(s, s).item()
            den = torch.dot(s, y).item()
            a_bb = num / den if den > 0 else 1.0  # update step size for next iteration

            g = g_new.clone()
            x = x_new.clone()

            err = torch.norm(self.project(x_new - g_new, *self.proj_args) - x_new).item()
            if err < tol:
                return x_new, True
            
            x, g = x_new, g_new  # advance iterate
        
        return x, False


    def _line_min(self, fun_jac, x0, f0, g0, p, f_b, args, amax=1, c1=0.0001, maxiter=10, verbal=False):
        """
        Perfoms line minimization to find optimal step size

        Args:
            fun_jac: objective function and gradient
            x0: current point
            f0: function value at x0
            g0: gradient value at x0
            p: search direction
            f_b: function value for armijio condition (typically f0)
            args: objective function arguments
            amax: maximum step size
            c1: Armijo condition parameter
            maxiter: maximum iterations
            verbal: whether to print information

        Returns:
            a: best step size
            f: objective function at best step size
            g: objective gradient at best step size
        """
        a_list = [0.0, float(amax)]
        f_list = [f0]
        g_list = [g0]
        i = 1
        while i <= maxiter:
            self.line_search_iteration = i
            if verbal:
                logging.log(f"     Line search iteration {i}: a = {a_list[i]:.6e}.", log_file=self.log_file)
            x_new = x0 + a_list[i] * p
            f_val, g_val = fun_jac(x_new, *args)
            if isinstance(f_val, torch.Tensor):
                f_val = f_val.item()
            f_list.append(f_val)
            g_list.append(g_val)

            if f_list[i] <= f_b + c1 * a_list[i] * torch.dot(g_list[0], p).item():
                return a_list[i], f_list[i], g_list[i]

            a_new = self._cubic_interp(
                a_list[i], f_list[i], torch.dot(g_list[i], p).item(),
                a_list[0], f_list[0], torch.dot(g_list[0], p).item()
            )
            if isinstance(a_new, complex) or a_new < 0:
                a_new = a_list[i] / 2.0 # just do back tracking if failure occurs
            a_list.append(a_new)
            i += 1

        return None, f_list[maxiter], g_list[maxiter]

    def _cubic_interp(self, an, fn, gn, ao, fo, go):
        """
        Obtain new step size by cubic interpolation

        Args:
            an: newest step size
            fn: fun evaluated at newest step size
            gn: grad evaluated at newest step size
            ao: previous step size
            fo: fun evaluated at previous step size
            go: grad evaluated at previous step size
        """
        # If error occurs, return step size an
        try:
            d1 = go + gn - 3.0 * (fo - fn) / (ao - an)
            d2 = (1.0 if (an - ao) >= 0 else -1.0) * ((d1 ** 2 - go * gn) ** 0.5)
            res = an - (an - ao) * ((gn + d2 - d1) / (gn - go + 2.0 * d2))
        except:
            res = an
        return res

    def _entropy(self, x):
        """
        Computes the Tsallis q-Entropy.
        """
        x_full = x.view(self.n_atoms, self.n_species)
        t = 0
        grad_t = torch.zeros_like(x_full)
        for i in range(self.n_atoms):
            if self.q==1:
                t += -torch.sum(x_full[i, :] * torch.log(x_full[i, :]+1e-10))
            else: 
                t += 1/(self.q-1) * (1-torch.sum(x_full[i, :]**self.q))
        return t