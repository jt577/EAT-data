###########################################################################
# Tsallis q entropy and its gradient
###########################################################################
import numpy as np

def Tsallis(w, S, q):
    """
    Computes the Tsallis q-Entropy

    Args:
        w: Probabilities (vectorized)
        S: Number of species
        q: q-factor for Entropy (q=1 is Shannon Entropy)

    Returns:
        res: The entropy
    """
    w_full = w.reshape(-1, S)
    N = w_full.shape[0]
    w_full = w.reshape(N, -1)
    res = 0
    for i in range(N):
        if q==1:
            res += -np.sum(w_full[i, :] * np.log(w_full[i, :]+1e-10))
        else: 
            res += 1/(q-1) * (1-np.sum(w_full[i, :]**q))

    return res

def grad_Tsallis(w, S, q):
    """
    Computes the gradient of the Tsallis q-Entropy

    Args:
        w: Probabilities (vectorized)
        S: Number of species
        q: q-factor for Entropy (q=1 is Shannon Entropy)

    Returns:
        res: Gradient of the entropy
    """
    w_full = w.reshape(-1, S)
    N = w_full.shape[0]
    grad_full = np.zeros(w_full.shape)
    for i in range(N):
        if q==1:
            grad_full[i, :] = -(np.log(w_full[i, :] + 1e-10) + 1)
        else: 
            grad_full[i, :] = -q/(q-1) * w_full[i, :]**(q-1)
    res = grad_full.flatten()
    return res