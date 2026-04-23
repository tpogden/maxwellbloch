# -*- coding: utf-8 -*-

"""Named time-function factories for use with the time-dependent solver.

Each function takes an integer index and returns a callable with signature
``f(t, args)`` where ``args`` is a dict of named parameters suffixed by the
index. For example, ``square(1)`` returns a function that reads
``{'on_1', 'off_1', 'ampl_1'}`` from ``args``.
"""

from collections.abc import Callable
from typing import Any

import numpy as np
from numpy import exp, log, pi, sqrt
from numpy import sinc as npsinc
from scipy.interpolate import interp1d

# Type alias for a time function returned by the factories below.
# The function takes a time value (scalar or array) and a parameter dict,
# and returns the corresponding amplitude (scalar or array).
TFunc = Callable[[float | np.ndarray, dict[str, Any]], float | np.ndarray]


def square(index: int) -> TFunc:
    """Return a square (top-hat) pulse time function.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``ampl_{index}``, ``on_{index}``,
        ``off_{index}`` from args. Returns ``ampl`` between ``on`` and
        ``off``, zero outside.
    """

    def func(t, args):
        on = args["on_" + str(index)]
        off = args["off_" + str(index)]
        ampl = args["ampl_" + str(index)]
        return ampl * (t >= on) * (t <= off)

    func.__name__ = "square_" + str(index)
    return func


def gaussian(index: int) -> TFunc:
    """Return a Gaussian pulse time function.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``fwhm_{index}``, ``centre_{index}``,
        and either ``ampl_{index}`` or ``n_pi_{index}`` from args. Raises
        ``KeyError`` if both or neither amplitude parameters are present.

    Notes:
        The amplitude can be set directly via ``ampl_{index}`` or indirectly
        via ``n_pi_{index}`` (desired pulse area in multiples of π).
    """

    def func(t, args):
        fwhm = args["fwhm_{}".format(index)]
        centre = args["centre_{}".format(index)]
        ampl_idx = "ampl_{}".format(index)
        n_pi_idx = "n_pi_{}".format(index)
        if ampl_idx in args:
            if n_pi_idx in args:
                raise KeyError("t_args can contain ampl or n_pi, not both.")
            else:
                ampl = args[ampl_idx]
        else:
            if n_pi_idx in args:
                n_pi = args[n_pi_idx]
                ampl = n_pi * sqrt(4.0 * pi * log(2) / fwhm**2) / (2 * pi)
            else:
                raise KeyError("t_args must contain ampl or n_pi.")
        return ampl * exp(-4 * log(2) * ((t - centre) / fwhm) ** 2)

    func.__name__ = "gaussian_{}".format(index)
    return func


def sech(index: int) -> TFunc:
    """Return a sech pulse time function.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``centre_{index}``, either
        ``ampl_{index}`` or ``n_pi_{index}``, and either ``width_{index}``
        or ``fwhm_{index}`` from args. Raises ``KeyError`` if conflicting
        or missing parameters are found.

    Notes:
        - The amplitude can be set directly via ``ampl_{index}`` or
          indirectly via ``n_pi_{index}`` (desired pulse area in multiples
          of π).
        - The width can be set directly via ``width_{index}`` or indirectly
          via ``fwhm_{index}`` (full-width at half maximum).
    """

    def sech_(t):
        return 2 / (exp(t) + exp(-t))

    SECH_FWHM_CONV = 1.0 / 2.6339157938

    def func(t, args):
        centre = args[f"centre_{index}"]
        width_idx = f"width_{index}"
        fwhm_idx = f"fwhm_{index}"
        ampl_idx = f"ampl_{index}"
        n_pi_idx = f"n_pi_{index}"
        if width_idx in args:
            if fwhm_idx in args:
                raise KeyError("t_args can contain width or fwhm, not both.")
            else:
                width = args[width_idx]
        else:
            if fwhm_idx in args:
                width = args[fwhm_idx] * SECH_FWHM_CONV
            else:
                raise KeyError("t_args must contain width or fwhm.")
        if ampl_idx in args:
            if n_pi_idx in args:
                raise KeyError("t_args can contain ampl or n_pi, not both.")
            else:
                ampl = args[ampl_idx]
        else:
            if n_pi_idx in args:
                n_pi = args[n_pi_idx]
                ampl = n_pi / width / (2 * pi)
            else:
                raise KeyError("t_args must contain ampl or n_pi.")
        return ampl * sech_((t - centre) / width)

    func.__name__ = "sech_" + str(index)
    return func


def ramp_on(index: int) -> TFunc:
    """Return a ramp-on time function.

    The pulse rises smoothly from zero using a half-Gaussian, then holds at
    ``ampl`` after the turn-on time ``on``.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``ampl_{index}``, ``fwhm_{index}``,
        and ``on_{index}`` from args.
    """

    def func(t, args):
        ampl = args["ampl_" + str(index)]
        fwhm = args["fwhm_" + str(index)]
        on = args["on_" + str(index)]
        return ampl * (exp(-4 * log(2) * ((t - on) / fwhm) ** 2) * (t <= on) + (t > on))

    func.__name__ = "ramp_on_" + str(index)
    return func


def ramp_off(index: int) -> TFunc:
    """Return a ramp-off time function.

    The pulse holds at ``ampl`` until the turn-off time ``off``, then falls
    smoothly to zero using a half-Gaussian.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``ampl_{index}``, ``fwhm_{index}``,
        and ``off_{index}`` from args.
    """

    def func(t, args):
        ampl = args["ampl_" + str(index)]
        fwhm = args["fwhm_" + str(index)]
        off = args["off_" + str(index)]
        return ampl * (
            exp(-4 * log(2) * ((t - off) / fwhm) ** 2) * (t >= off) + (t < off)
        )

    func.__name__ = "ramp_off_" + str(index)
    return func


def ramp_onoff(index: int) -> TFunc:
    """Return a ramp-on / ramp-off time function.

    The pulse rises smoothly, holds at ``ampl``, then falls smoothly. Built
    by composing :func:`ramp_on` and :func:`ramp_off`.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``ampl_{index}``, ``fwhm_{index}``,
        ``on_{index}``, and ``off_{index}`` from args.
    """
    _ramp_on = ramp_on(index)
    _ramp_off = ramp_off(index)

    def func(t, args):
        ampl = args["ampl_" + str(index)]
        return _ramp_on(t, args) + _ramp_off(t, args) - ampl

    func.__name__ = "ramp_onoff_" + str(index)
    return func


def ramp_offon(index: int) -> TFunc:
    """Return a ramp-off / ramp-on time function.

    The pulse starts at ``ampl``, dips smoothly to zero, then rises back to
    ``ampl``. Built by composing :func:`ramp_on` and :func:`ramp_off`.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``ampl_{index}``, ``fwhm_{index}``,
        ``off_{index}``, and ``on_{index}`` from args.
    """
    _ramp_on = ramp_on(index)
    _ramp_off = ramp_off(index)

    def func(t, args):
        return _ramp_on(t, args) + _ramp_off(t, args)

    func.__name__ = "ramp_offon_" + str(index)
    return func


def sinc(index: int) -> TFunc:
    """Return a sinc pulse time function.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``ampl_{index}`` and
        ``width_{index}`` from args.
    """

    def func(t, args):
        ampl = args["ampl_" + str(index)]
        width = args["width_" + str(index)]
        return ampl * npsinc(width * t) / sqrt(pi / 2.0)

    func.__name__ = "sinc_" + str(index)
    return func


def intp(index: int) -> TFunc:
    """Return an interpolated time function.

    Linearly interpolates a user-supplied ``(tlist, ylist)`` pair, returning
    zero outside the supplied range.

    Args:
        index: Integer suffix used to look up parameters in args.

    Returns:
        Callable ``f(t, args)`` reading ``tlist_{index}`` and
        ``ylist_{index}`` from args. Returns a Python complex scalar when
        called with scalar ``t`` (as required by QuTiP 5), or an array when
        called with array ``t``.
    """

    def func(t, args):
        tlist = args["tlist_" + str(index)]
        ylist = args["ylist_" + str(index)]
        yintp = interp1d(tlist, ylist, bounds_error=False, fill_value=0.0)
        result = yintp(t)
        # QuTiP 5 requires a Python scalar (not a 0-d ndarray) when called
        # with scalar t. When called with array t (e.g. init_Omegas_zt),
        # return the array as-is.
        return complex(result) if np.ndim(result) == 0 else result

    func.__name__ = "intp_" + str(index)
    return func
