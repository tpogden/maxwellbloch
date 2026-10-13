"""T0.5: finite-difference gradient harness."""

import jax
import jax.numpy as jnp
import numpy as np


def central_fd(f, x, idx, h):
    """Central FD of scalar f w.r.t. components ``idx`` of 1-D x, step h (absolute)."""
    x = np.asarray(x, dtype=float)
    out = np.empty(len(idx))
    for k, i in enumerate(idx):
        e = np.zeros_like(x)
        e[i] = h
        out[k] = (float(f(jnp.asarray(x + e))) - float(f(jnp.asarray(x - e)))) / (2 * h)
    return out


def fd_step_scan(f, x, idx, hs=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)):
    """FD at several steps; the best step is where successive estimates agree."""
    return {h: central_fd(f, x, idx, h) for h in hs}


def rel_err(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def check_gradient(loss, x, idx, hs=(1e-3, 1e-4, 1e-5)):
    """Compare jax.grad against central FD (best step) at components idx.
    Returns dict with AD values, FD per step, and rel error of the best step."""
    g = np.asarray(jax.jit(jax.grad(loss))(x))[idx]
    scan = fd_step_scan(jax.jit(loss), x, idx, hs)
    errs = {h: rel_err(g, v) for h, v in scan.items()}
    best = min(errs, key=errs.get)
    return {"ad": g, "fd": scan, "rel_err": errs, "best_h": best}
