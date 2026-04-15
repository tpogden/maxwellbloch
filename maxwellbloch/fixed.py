# -*- coding: utf-8 -*-

"""Functions for transforming solved results into the fixed (lab) frame of
reference, accounting for the finite speed of light."""

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def t_list(mb_solve, speed_of_light):
    """Return the time points shifted to the fixed (lab) frame of
    reference given a speed-of-light.

    Args:
        mb_solve: An MBSolve object
        speed_of_light: The speed of light in the system.

    Returns:
        Array of time values in the fixed frame of reference.
    """

    t_scale = 1.0 + mb_solve.z_max / (speed_of_light * mb_solve.t_max)
    return mb_solve.tlist * t_scale


def rabi_freq(mb_solve, field_idx, speed_of_light, part="real", interp_kind="linear"):
    """Return the field results shifted to the fixed (lab) frame of reference
    given a speed-of-light by interpolation. Can return the real part or
    abs value.

    Args:
        mb_solve: An MBSolve object
        field_idx: The field to return
        speed_of_light: The speed of light in the system
        part: Which part of the complex field ('real' or 'abs')
        interp_kind: The kind of spline interpolation to use ('linear',
            'cubic' or 'quintic')

    Returns:
        Array[num_fields, num_space_points, num_time_points] of field
        values in the fixed frame of reference.
    """

    if part == "abs":
        rabi_freq_zt = np.abs(mb_solve.Omegas_zt[field_idx])
    elif part == "real":
        rabi_freq_zt = np.real(mb_solve.Omegas_zt[field_idx])
    else:
        raise ValueError('Invalid part. Try "abs" or "real"')

    rabi_freq_intp = RegularGridInterpolator(
        (mb_solve.zlist, mb_solve.tlist),
        rabi_freq_zt,
        method=interp_kind,
        bounds_error=False,
        fill_value=0.0,
    )

    rabi_freq_fixed = np.zeros(mb_solve.Omegas_zt[field_idx].shape, dtype=float)

    for i, z_i in enumerate(mb_solve.zlist):
        t_shifted = t_list(mb_solve, speed_of_light) - z_i / speed_of_light
        pts = np.column_stack([np.full(len(t_shifted), z_i), t_shifted])
        rabi_freq_fixed[i] = rabi_freq_intp(pts)

    return rabi_freq_fixed


