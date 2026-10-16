"""Differentiable interpolation of sampled fields.

``PropagationSolver`` drives the Hamiltonian with ``qutip.coefficient(array,
tlist)``, which is a cubic *not-a-knot* spline (``scipy`` ``make_interp_spline``
with ``k=3, bc_type=None``). This module reproduces exactly that spline in JAX:
the knot derivatives are the solution of a fixed tridiagonal system whose
right-hand side is linear in the samples (O(n_t) in time and memory), and the
spline is then evaluated piecewise. Gradients flow into the samples.

JAX's own splines (natural, Hermite) are different interpolants and would not
match clerq.
"""

from __future__ import annotations

import jax.numpy as jnp
import numpy as np
from jax.lax.linalg import tridiagonal_solve
from scipy.interpolate import make_interp_spline


def not_a_knot_derivative_matrix(ts: np.ndarray) -> np.ndarray:
    """Return ``D`` (n x n) with ``D @ y`` the knot derivatives of the spline of ``y``.

    Dense, so O(n^2) in memory: this is the SciPy *reference* the banded
    :class:`CubicDrive` is tested against, not used at run time.
    """
    ts = np.asarray(ts, dtype=float)
    spline = make_interp_spline(ts, np.eye(len(ts)), k=3, bc_type=None)
    return spline.derivative()(ts)  # [i, j] = s_j'(t_i)


def _not_a_knot_tridiagonal(ts: np.ndarray):
    """The tridiagonal system ``A m = b`` for the knot derivatives ``m``.

    Returns ``(lower, diag, upper, first, last)``: the three bands of ``A``
    (``lower[i]`` is ``A[i+1, i]``, ``upper[i]`` is ``A[i, i+1]``), and the weights
    ``first = (w0, w1)`` and ``last = (w_-2, w_-1)`` such that ``b[0] = w0 * slope[0]
    + w1 * slope[1]`` and ``b[-1] = w_-2 * slope[-2] + w_-1 * slope[-1]``. The
    interior right-hand sides are ``3 (dx[i] slope[i-1] + dx[i-1] slope[i])``.
    This is the system SciPy's
    ``CubicSpline(bc_type="not-a-knot")`` solves, where the end conditions have
    been eliminated down to two terms so that the matrix stays tridiagonal.
    """
    dx = np.diff(ts)
    n = len(ts)
    diag = np.empty(n)
    diag[1:-1] = 2.0 * (dx[:-1] + dx[1:])
    diag[0], diag[-1] = dx[1], dx[-2]
    upper = np.empty(n - 1)
    upper[0], upper[1:] = ts[2] - ts[0], dx[:-1]
    lower = np.empty(n - 1)
    lower[:-1], lower[-1] = dx[1:], ts[-1] - ts[-3]
    d0, d1 = ts[2] - ts[0], ts[-1] - ts[-3]
    first = ((dx[0] + 2.0 * d0) * dx[1] / d0, dx[0] ** 2 / d0)
    last = (dx[-1] ** 2 / d1, (2.0 * d1 + dx[-1]) * dx[-2] / d1)
    return lower, diag, upper, first, last


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
        # O(n_t): the knot derivatives solve a fixed tridiagonal system.
        lower, diag, upper, self._first, self._last = _not_a_knot_tridiagonal(ts_np)
        # jax's tridiagonal_solve wants dl, du padded to the length of d
        self._dl = jnp.asarray(np.concatenate([[0.0], lower]))
        self._d = jnp.asarray(diag)
        self._du = jnp.asarray(np.concatenate([upper, [0.0]]))

    def knot_derivatives(self, y):
        """``s'(t_i)`` of the spline through ``y`` (shape ``(..., n_t)``), same shape."""
        y = jnp.asarray(y)
        slope = jnp.diff(y, axis=-1) / self.h
        mid = 3.0 * (self.h[1:] * slope[..., :-1] + self.h[:-1] * slope[..., 1:])
        b0 = self._first[0] * slope[..., 0] + self._first[1] * slope[..., 1]
        b1 = self._last[0] * slope[..., -2] + self._last[1] * slope[..., -1]
        rhs = jnp.concatenate([b0[..., None], mid, b1[..., None]], axis=-1)
        batch = rhs.shape[:-1]
        cols = rhs.reshape(-1, rhs.shape[-1]).T  # (n_t, number of right-hand sides)
        m = tridiagonal_solve(self._dl, self._d, self._du, cols)
        return m.T.reshape(*batch, -1)

    def coefficients(self, y):
        """Per-interval ``(a, b, c, d)`` of ``d + c s + b s^2 + a s^3``, ``s = t - t_i``."""
        y = jnp.asarray(y)
        m = self.knot_derivatives(y)
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
