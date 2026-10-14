"""Differentiable interpolation of sampled fields.

``PropagationSolver`` drives the Hamiltonian with ``qutip.coefficient(array,
tlist)``, which is a cubic *not-a-knot* spline (``scipy`` ``make_interp_spline``
with ``k=3, bc_type=None``). This module reproduces exactly that spline in JAX:
the knot derivatives are a fixed linear map ``D @ y`` of the samples (``D`` is
computed once with SciPy), and the spline is then evaluated piecewise. Gradients
flow into the samples.

JAX's own splines (natural, Hermite) are different interpolants and would not
match clerq.
"""

from __future__ import annotations

import jax.numpy as jnp
import numpy as np
from scipy.interpolate import make_interp_spline


def not_a_knot_derivative_matrix(ts: np.ndarray) -> np.ndarray:
    """Return ``D`` (n x n) with ``D @ y`` the knot derivatives of the spline of ``y``."""
    ts = np.asarray(ts, dtype=float)
    spline = make_interp_spline(ts, np.eye(len(ts)), k=3, bc_type=None)
    return spline.derivative()(ts)  # [i, j] = s_j'(t_i)


class CubicDrive:
    """Cubic not-a-knot interpolation of samples on a fixed time grid.

    Calling an instance with samples ``y`` of shape ``(..., n_t)`` returns a
    function ``f(t)`` that evaluates the splines at a scalar time, with shape
    ``(...,)``. Leading axes (for example one per field) are interpolated
    together.

    Args:
        ts: strictly increasing sample times, at least four.

    Raises:
        ValueError: if the times are not strictly increasing or there are fewer
            than four.
    """

    def __init__(self, ts):
        ts_np = np.asarray(ts, dtype=float)
        if ts_np.ndim != 1 or len(ts_np) < 4:
            raise ValueError("A cubic spline needs a 1-D grid of at least 4 times.")
        if np.any(np.diff(ts_np) <= 0.0):
            raise ValueError("Sample times must be strictly increasing.")
        self.ts = jnp.asarray(ts_np)
        self.h = jnp.diff(self.ts)
        # Dense n_t x n_t; a banded solve would scale better for very long grids.
        self.D = jnp.asarray(not_a_knot_derivative_matrix(ts_np))

    def coefficients(self, y):
        """Per-interval ``(a, b, c, d)`` of ``d + c s + b s^2 + a s^3``, ``s = t - t_i``."""
        y = jnp.asarray(y)
        m = jnp.einsum("ij,...j->...i", self.D, y)
        slope = jnp.diff(y, axis=-1) / self.h
        a = (m[..., :-1] + m[..., 1:] - 2.0 * slope) / self.h**2
        b = (3.0 * slope - 2.0 * m[..., :-1] - m[..., 1:]) / self.h
        return a, b, m[..., :-1], y[..., :-1]

    def __call__(self, y):
        a, b, c, d = self.coefficients(y)
        ts = self.ts
        last = len(ts) - 2

        def f(t):
            i = jnp.clip(jnp.searchsorted(ts, t, side="right") - 1, 0, last)
            s = t - ts[i]
            pick = lambda arr: jnp.take(arr, i, axis=-1)  # noqa: E731
            return pick(d) + s * (pick(c) + s * (pick(b) + s * pick(a)))

        return f
