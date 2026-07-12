import math

import torch


def softplus(E, beta):
    """
    Softplus cost function plus derivative factor.

    Supports either Python scalars or torch scalars.
    """
    if torch.is_tensor(E):
        beta_t = torch.as_tensor(beta, dtype=E.dtype, device=E.device)
        y = beta_t * E
        f = torch.nn.functional.softplus(y) / beta_t
        g = torch.sigmoid(y)
        return f, g

    y = float(beta) * float(E)
    f = math.log1p(math.exp(y)) / float(beta)
    g = math.exp(y) / (1.0 + math.exp(y))
    return f, g


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
        grad_ent = -Tsallis_q * X ** (Tsallis_q - 1.0) / (Tsallis_q - 1.0)

    return ent, grad_ent


def max_tsallis_entropy(num_states, Tsallis_q):
    if num_states <= 1:
        return 0.0
    if Tsallis_q == 1.0:
        return math.log(float(num_states))
    return (1.0 - float(num_states) ** (1.0 - Tsallis_q)) / (Tsallis_q - 1.0)


def entropy_ceiling(n_sites, n_elements, Tsallis_q, n_null_sites=0):
    species_max = max_tsallis_entropy(n_elements, Tsallis_q)
    null_max = max_tsallis_entropy(n_elements + 1, Tsallis_q)
    n_null_sites = int(n_null_sites)
    n_normal_sites = int(n_sites) - n_null_sites
    return n_normal_sites * species_max + n_null_sites * null_max


def null_total_entropy(W, Tsallis_q, null_mask, eps=1e-10):
    """
    Total entropy for a mix of normal EAT rows and null rows.

    Normal rows contribute the usual categorical entropy.
    Null rows contribute one categorical entropy over
      [p_species_1, ..., p_species_n, p_empty]
    where p_empty = 1 - sum_i p_species_i.
    """
    device = W.device
    null_mask = torch.as_tensor(null_mask, device=device, dtype=torch.bool).view(-1)
    if null_mask.numel() != W.shape[0]:
        raise ValueError("null_mask must have one entry per row in W")

    total_species_entropy = torch.zeros((), dtype=W.dtype, device=device)
    total_null_entropy = torch.zeros((), dtype=W.dtype, device=device)

    normal_mask = ~null_mask
    if torch.any(normal_mask):
        total_species_entropy = total_species_entropy + entropy(W[normal_mask], Tsallis_q, eps=eps)[0]

    if torch.any(null_mask):
        Wg = W[null_mask]
        empty = torch.clamp(1.0 - Wg.sum(dim=1, keepdim=True), min=0.0)
        null_probs = torch.cat([Wg, empty], dim=1)
        total_null_entropy = total_null_entropy + entropy(null_probs, Tsallis_q, eps=eps)[0]

    total_entropy = total_species_entropy + total_null_entropy
    return total_entropy, total_species_entropy, total_null_entropy
