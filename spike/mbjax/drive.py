"""Phase 2: differentiable drive from a sampled array.

clerq interpolates the propagating field with QuTiP's ``coefficient(array, tlist)``,
i.e. a cubic spline from ``scipy.interpolate.make_interp_spline(k=3, bc_type=None)``
(not-a-knot). We reproduce exactly that spline: the knot derivatives are a *fixed linear
map* of the samples, ``m = D @ y`` (D is computed once with scipy), and the spline is then
evaluated piecewise from (y, m). Everything downstream is JAX, so gradients flow into y.

``LinearDrive`` (np.interp) is provided for comparison only.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
from scipy.interpolate import make_interp_spline

jax.config.update("jax_enable_x64", True)


def not_a_knot_derivative_matrix(ts: np.ndarray) -> np.ndarray:
    """D (n x n) with ``D @ y`` = knot first-derivatives of the not-a-knot cubic spline."""
    ts = np.asarray(ts, dtype=float)
    spl = make_interp_spline(ts, np.eye(len(ts)), k=3, bc_type=None)
    return spl.derivative()(ts)  # [i, j] = s_j'(t_i); m_i = sum_j D[i,j] y_j


class CubicDrive:
    """Callable ``t -> value`` for samples ``y`` on knots ``ts`` (cubic, not-a-knot)."""

    def __init__(self, ts, D=None):
        self.ts = jnp.asarray(ts)
        self.D = jnp.asarray(not_a_knot_derivative_matrix(ts) if D is None else D)
        self.h = jnp.diff(self.ts)

    def coefficients(self, y):
        """Per-interval polynomial coefficients (a, b, c, d): d + c s + b s^2 + a s^3,
        s = t - ts[i]."""
        m = self.D @ y
        dy = jnp.diff(y)
        h = self.h
        c = m[:-1]
        b = (3.0 * dy / h - 2.0 * m[:-1] - m[1:]) / h
        a = (m[:-1] + m[1:] - 2.0 * dy / h) / h**2
        return a, b, c, y[:-1]

    def __call__(self, y):
        """Return an evaluator ``f(t)`` for samples y."""
        a, b, c, d = self.coefficients(y)
        ts = self.ts

        def f(t):
            i = jnp.clip(jnp.searchsorted(ts, t, side="right") - 1, 0, len(ts) - 2)
            s = t - ts[i]
            return d[i] + s * (c[i] + s * (b[i] + s * a[i]))

        return f


def linear_drive(ts):
    ts = jnp.asarray(ts)
    return lambda y: (lambda t: jnp.interp(t, ts, y))
