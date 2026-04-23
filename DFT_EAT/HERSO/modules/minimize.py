###################################################################################################
# Functions for Projected Quasi Newton Minimization
###################################################################################################

# imports
import numpy as np
import modules.bsruncalc as bs


def create_diagonal_vectors(N, S):
    I = np.eye(N)
    result = np.repeat(I, S, axis=0)
    return result

def create_block_matrix(H, I):
    zero_matrix = np.zeros((I.shape[1], I.shape[1]))
    big_matrix = np.block([[H, I], [I.T, zero_matrix]])
    return big_matrix

def solve_linear_system(M, Hx0, g, N):
    one_vec = np.ones(N)
    r = (-g + Hx0).flatten()
    c = np.concatenate((r, one_vec))
    z = np.linalg.solve(M, c)
    x = z[:-N]
    lambd = z[-N:]
    return x, lambd

def solve_quadratic_form(x0, g, H, N, S, details_file_path, max_iter=10000):
    with open(details_file_path, 'a') as file:
        file.write(f"Optimizing BFGS quadratic form subproblem...\n")
    H = 0.5 * (H + H.T)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    eigenvalues[eigenvalues <= 0] = 1e-10
    H = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T

    I = create_diagonal_vectors(N, S)
    Hx0 = H @ x0
    A = np.arange(N * S)
    numit = 1
    converged = False
    tolerance = 1e-6  # Adjusted tolerance

    while not converged and numit <= max_iter:
        x0_A = x0[A]
        g_A = g[A]
        H_AA = H[np.ix_(A, A)]
        I_A = I[A, :]
        Hx0_A = Hx0[A]

        M = create_block_matrix(H_AA, I_A)
        x_A, lambd = solve_linear_system(M, Hx0_A, g_A, N)

        x = np.zeros(N * S)
        x[A] = x_A

        grad_L = H @ (x - x0) + g + I @ lambd
        mu = grad_L
        mu[A] = 0

        not_in_A = np.setdiff1d(np.arange(N * S), A)
        mu_not_in_A = mu[not_in_A]

        indices_to_add = not_in_A[mu_not_in_A < -tolerance]
        x_A = x[A]
        indices_to_remove = A[x_A < -tolerance]

        if len(indices_to_add) == 0 and len(indices_to_remove) == 0:
            converged = True
        else:
            if len(indices_to_remove) > 0:
                min_x_index = indices_to_remove[np.argmin(x_A[x_A < -tolerance])]
                indices_to_remove = [min_x_index]
                with open(details_file_path, 'a') as file:
                    file.write(f"    Iteration {numit}: Removing variable {min_x_index} from active set.\n")
                x[min_x_index] = 0
                A = np.setdiff1d(A, indices_to_remove)
            if len(indices_to_add) > 0:
                min_mu_index = indices_to_add[np.argmin(mu_not_in_A[mu_not_in_A < -tolerance])]
                indices_to_add = [min_mu_index]
                with open(details_file_path, 'a') as file:
                    file.write(f"    Iteration {numit}: Adding variable {min_mu_index} to active set.\n")
                A = np.union1d(A, indices_to_add)
            numit += 1

    if not converged:
        with open(details_file_path, 'a') as file:
            file.write("    Warning: Maximum iterations reached without convergence.\n")

    value = 0.5 * (x - x0).T @ H @ (x - x0) + g.T @ (x - x0)
    with open(details_file_path, 'a') as file:
        file.write(f'    Successfully optimized after {numit - 1} iterations.\n    Df (must be < 0) = {value}\n')

    return x


def project_to_simplex_fast(x, N):
    phi = x.reshape(N, -1)
    # Sort phi in descending order
    phi_sorted = np.sort(phi, axis=-1)[:, ::-1]
    
    # Calculate the cumulative sum of sorted phi
    cumsum_phi = np.cumsum(phi_sorted, axis=-1)
    
    # Determine the number of elements
    n = phi.shape[-1]
    
    # Find the threshold by iterating over each grid point
    rho = np.arange(1, n+1)
    rho_matrix = np.tile(rho, (phi.shape[0], 1))
    condition = phi_sorted - (cumsum_phi - 1) / rho_matrix > 0
    
    # Find the largest k where the condition is True
    k = np.max(np.where(condition, rho_matrix, 0), axis=-1)
    
    # Compute the threshold theta
    theta = (cumsum_phi[np.arange(phi.shape[0]), k-1] - 1) / k
    
    # Project onto the simplex
    phi_projected = np.maximum(phi - theta[:, np.newaxis], 0)
    
    return phi_projected.reshape(-1)

def PBFGS(fun, x0, jac, args, S,  folder, progress_file_path, details_file_path, callback, B=None, project=project_to_simplex_fast, maxit=10, tol=1e-3):
    """
    Perform projected BFGS

    Args:
        fun: Function to minimize
        x0: Initial guess
        jac: Gradient of function
        args: Extra args passed to fun and jac
        S: Number of species
        folder: Folder from which we extract energies and gradients
        progress_file_path: Path to progress file
        details_file_path: Path to details file
        callback: Callback function
        B: Initial quasi-Hessian
        project: Projection back onto nearest point on valid domain
        maxit: Maximum iterations
        tol: Tolerance of convergence

    Returns:
        x: Final solution
        B: Final Hessian
        unique_folder: Folder from which we extract energies and gradients
    """
    unique_folder = folder
    x_full = x0.reshape(-1, S)
    N = x_full.shape[0]
    x = x0
    I = np.eye(len(x))
    f_new = None
    g_new = None
    k = 1
    while k <= maxit:
        if f_new is not None and g_new is not None:
            f = f_new
            g = g_new
        else:
            # Perform JDFTx calculation
            unique_folder = bs.perform_calc(x, *args)
            f = fun(x, unique_folder, *args)
            g = jac(x, unique_folder, *args)
        # Write to min_progress file
        callback(x, unique_folder)
        # Check for convergence
        if np.linalg.norm(project(x-g, N) - x) < tol:
            with open(details_file_path, 'a') as file:
                file.write(f"\nPBFGS converged after {k} iterations: fun = {f}, tol = {np.linalg.norm(project(x-g, N) - x)}\nx = {x}\n")
            break
        # Print current iteration summary
        with open(details_file_path, 'a') as file:
            file.write(f"\nPBFGS iteration {k}: fun = {f}, tol = {np.linalg.norm(project(x-g, N) - x)}\n")
        if B is None:
            B = 5 *  np.linalg.norm(g) * I # Initial quasi-Hessian
        # Obtain line search direction
        x_star = solve_quadratic_form(x0=x, g=g, H=B, N=N, S=S, details_file_path=details_file_path)
        x_star = project(x_star, N) # project result to be safe
        d = x_star - x
        # Choose step by line search
        alpha = 1
        alpha, f_new, g_new, unique_folder = line_min(fun=fun, jac=jac, x0=x, f0=f, g0=g, p=d, args=args, folder=unique_folder, details_file_path=details_file_path, amax=alpha, c1=0.0, maxiter=4)
        
        # If line search failed, reset Hessian
        if alpha is None:
            with open(details_file_path, 'a') as file:
                file.write(f"Line search failed, presumably because of bad curvature. Increasing entropy parameter.\n")
            return x, B, unique_folder
        
        # Update parameters
        x += alpha * d
        x = project(x, N) # project again to make sure no negative values
        s = alpha * d
        y = g_new - g

        # Damped BFGS Update
        if np.dot(s, y) >= 0.2 * s @ B @ s:
            theta = 1
        else:
            theta = 0.8 * s @ B @ s / (s @ B @ s - np.dot(s, y))
        r = theta * y + (1-theta) * B @ s
        B = B - B @ np.outer(s, s) @ B / (s @ B @ s) + np.outer(r, r) / np.dot(s, r) 
        # Write result to file if at maximum iterations
        if k == maxit:
            with open(details_file_path, 'a') as file:
                file.write(f"\nPBFGS did not converge after {k} iterations: fun = {f_new}, tol = {np.linalg.norm(project(x-g, N) - x)}\n")
        k += 1
    return x, B, unique_folder

def line_min(fun, jac, x0, f0, g0, p, args, folder, details_file_path, amax=1, c1=0.0001, maxiter=10):
    """
    Perfoms line minimization to find optimal step size

    Args:
        fun: objective function
        jac: objective gradient
        x0: current point
        f0: function value at x0
        g0: gradient value at x0
        p: search direction
        args: objective function arguments
        folder: Folder from which we extract energies and gradients
        details_file_path: Path to details file
        amax: maximum step size
        c1: Armijo condition parameter
        maxiter: maximum iterations

    Returns:
        a: best step size
        f: objective function at best step size
        g: objective gradient at best step size
        unique_folder: Folder from which we extract energies and gradients

    """
    unique_folder = folder
    a = []
    a.append(0) # first value of alpha is 0
    a.append(amax) # second value is amax
    f = []
    f.append(f0)
    g = []
    g.append(g0)
    i = 1
    with open(details_file_path, 'a') as file:
        file.write(f'Performing line minimization...\n')
    while i<=maxiter:
        # Perform JDFTx calculation
        with open(details_file_path, 'a') as file:
            file.write(f'   ->Iteration {i}\n')
        unique_folder = bs.perform_calc(x0 + a[i] * p, *args)      
        f.append(fun(x0 + a[i] * p, unique_folder, *args))
        g.append(jac(x0 + a[i] * p, unique_folder, *args))
        if f[i] <= f[0] + c1 * a[i] * np.dot(g[0], p):
            with open(details_file_path, 'a') as file:
                file.write(f'Line minimization succeded after {i} iterations.\n')
            return a[i], f[i], g[i], unique_folder
        else:
            a.append(cubic_interp(a[i], f[i], np.dot(g[i], p), a[0], f[0], np.dot(g[0], p)))
        i += 1
    with open(details_file_path, 'a') as file:
        file.write(f'Line minimization failed to find minimizing step after {i} iterations\n')
        return None, f[maxiter], g[maxiter], unique_folder

def cubic_interp(an, fn, gn, ao, fo, go):
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
    d1 = go + gn - 3*(fo - fn)/(ao - an)
    d2 = np.sign(an - ao) * np.sqrt(d1**2 - go * gn)
    res = an - (an-ao) * ((gn + d2 - d1)/(gn - go + 2*d2))
    return res
