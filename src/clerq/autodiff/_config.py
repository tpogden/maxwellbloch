"""Import and precision helpers for the optional autodiff backend."""

import importlib

_INSTALL_HINT = (
    "clerq.autodiff needs the optional JAX dependencies: pip install 'clerq[jax]'"
)


def require_jax() -> None:
    """Raise an ImportError with an install hint if the JAX stack is missing.

    Raises:
        ImportError: if ``jax`` or ``diffrax`` cannot be imported.
    """
    for name in ("jax", "diffrax"):
        try:
            importlib.import_module(name)
        except ImportError as e:
            raise ImportError(f"{_INSTALL_HINT} (could not import {name!r})") from e


def float_dtype():
    """The floating-point dtype JAX computes in: float64 with ``jax_enable_x64`` on, else float32."""
    import jax.numpy as jnp

    return jnp.result_type(float)


def complex_dtype():
    """The complex dtype matching :func:`float_dtype` (complex128 or complex64)."""
    import jax.numpy as jnp

    return jnp.result_type(complex)


def default_tolerances(
    rtol: float | None = None, atol: float | None = None
) -> tuple[float, float]:
    """Step-size tolerances of the adaptive solver, defaulting by precision.

    Float64 uses ``(1e-8, 1e-10)``. Float32 cannot resolve those (its machine epsilon is
    1.2e-7), and a controller asked for less than it can deliver takes many extra steps, so
    it defaults to ``(1e-6, 1e-9)``, which reproduces float64 results to about 1e-6 on the
    benchmark problems. An explicit value always wins.
    """
    double = float_dtype().itemsize == 8
    if rtol is None:
        rtol = 1e-8 if double else 1e-6
    if atol is None:
        atol = 1e-10 if double else 1e-9
    return float(rtol), float(atol)


def require_float64() -> None:
    """Raise if JAX is not computing in 64-bit precision.

    ``clerq.autodiff`` itself no longer requires this: it computes in whichever precision
    JAX is configured for (float32 by default, float64 with ``jax_enable_x64``). Call it
    where float64 is needed, for example for finite-difference reference checks.

    The global JAX configuration is never changed on the user's behalf: it
    would silently alter every other JAX computation in the process.

    Raises:
        RuntimeError: if ``jax_enable_x64`` is off.
    """
    import jax.numpy as jnp

    if jnp.asarray(0.0).dtype != jnp.float64:
        raise RuntimeError(
            "clerq.autodiff requires 64-bit precision, which JAX disables by "
            "default. Enable it before creating any arrays: "
            "jax.config.update('jax_enable_x64', True), or set the environment "
            "variable JAX_ENABLE_X64=1."
        )
