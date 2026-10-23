"""Helpers for parametrising free-form control fields for gradient-based design.

A control field given as one free value per time sample is a bad parameter set:
the optimiser is free to use features at the scale of the time grid (near the Nyquist
frequency), where the result depends on the grid and the interpolation rather than on
the physics, and can report efficiencies no physical field could reach. Expanding the
control in a few smooth functions removes that freedom, and also gives far fewer
parameters.

The pieces compose with ``AutodiffPropagation``: pass ``basis @ coefficients`` (bounded
if needed) as a row of ``omega0``, and differentiate with respect to the coefficients::

    t = ad.tlist
    B = controls.gaussian_basis(t, 12)
    a0 = controls.fit_coefficients(B, controls.unbounded(initial, 0.0, 5.0))

    def field(a):                       # cyclic Rabi frequency samples in [0, 5]
        return controls.bounded(B @ a, 0.0, 5.0)
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
from numpy.typing import ArrayLike


def gaussian_basis(
    t: ArrayLike,
    n: int,
    lo: float | None = None,
    hi: float | None = None,
    width: float | None = None,
) -> np.ndarray:
    """Gaussian bumps on a regular grid of centres, evaluated at the times ``t``.

    Args:
        t: sample times, shape ``(n_t,)``.
        n: number of basis functions (at least 2).
        lo: first centre; default the first time.
        hi: last centre; default the last time. Use a ``hi`` below the last time to
            confine a control to an earlier window.
        width: standard deviation of each bump; default the centre spacing, which
            makes neighbouring bumps overlap enough to be smooth.

    Returns:
        ``B``, shape ``(n_t, n)``, so that a control is ``B @ coefficients``.
    """
    t = np.asarray(t, dtype=float)
    if n < 2:
        raise ValueError(f"n must be at least 2, got {n}.")
    lo = float(t[0]) if lo is None else float(lo)
    hi = float(t[-1]) if hi is None else float(hi)
    if not hi > lo:
        raise ValueError(f"hi ({hi}) must be greater than lo ({lo}).")
    centres = np.linspace(lo, hi, n)
    sigma = (centres[1] - centres[0]) if width is None else float(width)
    return np.exp(-0.5 * ((t[:, None] - centres[None, :]) / sigma) ** 2)


def fit_coefficients(basis: ArrayLike, target: ArrayLike) -> np.ndarray:
    """Least-squares coefficients ``a`` with ``basis @ a`` closest to ``target``.

    Use it to start an optimisation from a known waveform. The fit is only as good as
    the basis can represent (a sharp edge is not represented by wide bumps).
    """
    return np.linalg.lstsq(
        np.asarray(basis), np.asarray(target, dtype=float), rcond=None
    )[0]


def bounded(x: ArrayLike, lower: float, upper: float) -> jax.Array:
    """Map real values to ``(lower, upper)`` with a sigmoid (JAX-differentiable).

    Keeps an amplitude within a budget (for example a peak Rabi frequency) without
    clipping, which would have zero gradient outside the range.
    """
    return lower + (upper - lower) * jax.nn.sigmoid(jnp.asarray(x))


def unbounded(
    y: ArrayLike, lower: float, upper: float, margin: float = 0.02
) -> np.ndarray:
    """Inverse of :func:`bounded`; values are first clipped ``margin`` inside the range.

    The clip keeps the logit finite and away from the flat tails of the sigmoid,
    where the gradient is almost zero.
    """
    if not 0.0 <= margin < 0.5:
        raise ValueError(f"margin must be in [0, 0.5), got {margin}.")
    frac = (np.asarray(y, dtype=float) - lower) / (upper - lower)
    frac = np.clip(frac, margin, 1.0 - margin)
    return np.log(frac / (1.0 - frac))
