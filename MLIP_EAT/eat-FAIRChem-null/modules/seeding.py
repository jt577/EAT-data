import random
import torch
import numpy as np

# -----------------------------
# RNG state helpers
# -----------------------------
def save_rng_state():
    st = {
        "py": random.getstate(),
        "np": np.random.get_state(),
        "torch_cpu": torch.get_rng_state(),
    }
    if torch.cuda.is_available():
        st["torch_cuda"] = torch.cuda.get_rng_state_all()
    return st


def load_rng_state(st):
    random.setstate(st["py"])
    np.random.set_state(st["np"])
    torch.set_rng_state(st["torch_cpu"])
    if torch.cuda.is_available() and "torch_cuda" in st:
        torch.cuda.set_rng_state_all(st["torch_cuda"])