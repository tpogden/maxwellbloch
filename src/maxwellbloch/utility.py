import numpy as np
from scipy import interpolate


def maxwell_boltzmann(v: np.ndarray, fwhm: float, offset: float = 0.0) -> np.ndarray:
    """Maxwell-Boltzmann probability distribution.

    Args:
        v: velocity (or detuning) values.
        fwhm: full-width at half-maximum of the distribution.
        offset: centre of the distribution (default 0).

    Returns:
        Probability density at each value in ``v``.
    """
    return 1.0 / (fwhm * np.sqrt(np.pi)) * np.exp(-(((v - offset) / fwhm) ** 2))


def half_max_roots(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Return the half-maximum value and the two roots where ``y`` crosses it.

    Args:
        x: 1-D array of independent variable values (must be monotonically
            increasing and have enough points for spline fitting).
        y: 1-D array of dependent variable values, same length as ``x``.

    Returns:
        Tuple ``(half_max, r1, r2)`` where ``half_max = max(y) / 2`` and
        ``r1 < r2`` are the two x-values where ``y == half_max``.
    """
    half_max = np.max(y) / 2
    spline = interpolate.UnivariateSpline(x, y - half_max, s=0)
    r1, r2 = spline.roots()
    return half_max, r1, r2


def full_width_at_half_max(x: np.ndarray, y: np.ndarray) -> float:
    """Return the full width at half maximum (FWHM) of a peak.

    Args:
        x: 1-D array of independent variable values (must be monotonically
            increasing and have enough points for spline fitting).
        y: 1-D array of dependent variable values, same length as ``x``.

    Returns:
        The FWHM as a scalar: the distance between the two points where
        ``y`` equals half its maximum value.
    """
    half_max, r1, r2 = half_max_roots(x, y)
    return r2 - r1
