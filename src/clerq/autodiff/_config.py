"""Import and precision guards for the optional autodiff backend."""

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


def require_float64() -> None:
    """Raise if JAX is not computing in 64-bit precision.

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
