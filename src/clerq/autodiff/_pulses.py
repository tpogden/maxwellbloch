"""JAX-traceable ports of clerq's pulse functions (``clerq.t_funcs``).

``PropagationSolver`` drives the first slice with the analytic ``rabi_freq_t_func``
of each field and later slices with an interpolation of the propagated field. The
analytic pulse has to be traceable for JAX, so the functions that clerq ships are
re-implemented here with identical formulas. Parameter names are those of the
JSON problem definitions (``ampl``, ``centre``, ``fwhm``, ``width``, ``n_pi``,
``on``, ``off``); clerq stores them with an ``_{index}`` suffix, which is stripped.

Unsupported pulse functions raise ``NotImplementedError``.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import jax.numpy as jnp

SECH_FWHM_CONV = 1.0 / 2.6339157938

SUPPORTED = (
    "gaussian",
    "sech",
    "square",
    "ramp_on",
    "ramp_off",
    "ramp_onoff",
    "ramp_offon",
)


def pulse_name(t_func: Callable) -> str:
    """Name of a clerq pulse function, from its ``__name__`` (``"gaussian_0"``)."""
    return t_func.__name__.rsplit("_", 1)[0]


def strip_index(args: dict) -> dict:
    """Remove clerq's ``_{index}`` suffix from the keys of a t_args dict."""
    return {k.rsplit("_", 1)[0]: v for k, v in args.items()}


def _need(args: dict, name: str, *keys: str):
    missing = [k for k in keys if k not in args]
    if missing:
        raise ValueError(f"Pulse {name!r} needs the argument(s) {missing}.")


def _gaussian_amplitude(args: dict):
    if "ampl" in args:
        if "n_pi" in args:
            raise KeyError("t_args can contain ampl or n_pi, not both.")
        return args["ampl"]
    _need(args, "gaussian", "n_pi")
    return (
        args["n_pi"]
        * math.sqrt(4 * math.pi * math.log(2) / args["fwhm"] ** 2)
        / (2 * math.pi)
    )


def pulse_function(name: str, args: dict) -> Callable:
    """Return ``f(t)``, the pulse ``name`` with (suffix-free) ``args``.

    Args:
        name: one of :data:`SUPPORTED`.
        args: pulse parameters; values may be JAX tracers.

    Returns:
        A function of a scalar or array time, built from JAX operations.

    Raises:
        NotImplementedError: if the pulse function is not supported.
        ValueError: if a required parameter is missing.
    """
    k4 = 4.0 * math.log(2.0)

    if name == "gaussian":
        _need(args, name, "fwhm", "centre")
        ampl = _gaussian_amplitude(args)
        return lambda t: (
            ampl * jnp.exp(-k4 * ((t - args["centre"]) / args["fwhm"]) ** 2)
        )

    if name == "sech":
        _need(args, name, "centre")
        if "width" in args:
            if "fwhm" in args:
                raise KeyError("t_args can contain width or fwhm, not both.")
            width = args["width"]
        else:
            _need(args, name, "fwhm")
            width = args["fwhm"] * SECH_FWHM_CONV
        if "ampl" in args:
            ampl = args["ampl"]
        else:
            _need(args, name, "n_pi")
            ampl = args["n_pi"] / width / (2.0 * math.pi)
        centre = args["centre"]

        def sech(t):
            x = (t - centre) / width
            return ampl * 2.0 / (jnp.exp(x) + jnp.exp(-x))

        return sech

    if name == "square":
        _need(args, name, "ampl", "on", "off")
        return lambda t: (
            args["ampl"] * jnp.where((t >= args["on"]) & (t <= args["off"]), 1.0, 0.0)
        )

    if name in ("ramp_on", "ramp_off", "ramp_onoff", "ramp_offon"):
        need = {
            "ramp_on": ("ampl", "fwhm", "on"),
            "ramp_off": ("ampl", "fwhm", "off"),
            "ramp_onoff": ("ampl", "fwhm", "on", "off"),
            "ramp_offon": ("ampl", "fwhm", "on", "off"),
        }[name]
        _need(args, name, *need)
        ampl, fwhm = args["ampl"], args["fwhm"]

        def ramp_on(t):
            on = args["on"]
            return ampl * (jnp.exp(-k4 * ((t - on) / fwhm) ** 2) * (t <= on) + (t > on))

        def ramp_off(t):
            off = args["off"]
            return ampl * (
                jnp.exp(-k4 * ((t - off) / fwhm) ** 2) * (t >= off) + (t < off)
            )

        if name == "ramp_on":
            return ramp_on
        if name == "ramp_off":
            return ramp_off
        if name == "ramp_onoff":
            return lambda t: ramp_on(t) + ramp_off(t) - ampl
        return lambda t: ramp_on(t) + ramp_off(t)

    raise NotImplementedError(
        f"Pulse function {name!r} is not supported by clerq.autodiff; "
        f"supported: {', '.join(SUPPORTED)}."
    )


def jump_times(name: str, args: dict) -> tuple:
    """Times at which the pulse is discontinuous, for step-size control."""
    if name == "square":
        return (args["on"], args["off"])
    return ()
