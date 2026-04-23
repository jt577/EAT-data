import numpy as np
import torch
import math

def softplus(E, beta):
    """
    Softplus formation energy cost function and gradient factor.
    """
    f = math.log(1 + math.exp(beta * E)) / beta
    g = math.exp(beta * E) / (1 + math.exp(beta * E))
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
        grad_ent = -Tsallis_q * X**(Tsallis_q - 1.0) / (Tsallis_q - 1.0)

    return ent, grad_ent