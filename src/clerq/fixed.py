# -*- coding: utf-8 -*-

"""Functions for transforming solved results into the fixed (lab) frame of
reference, accounting for the finite speed of light."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from clerq._deprecation import renamed_arg

if TYPE_CHECKING:
    from clerq.propagation import PropagationSolver


@renamed_arg("mb_solve", "solver")
def t_list(solver: PropagationSolver, speed_of_light: float) -> np.ndarray:
    """Return the time points shifted to the fixed (lab) frame of
    reference given a speed-of-light.

    Args:
        solver: A PropagationSolver object
        speed_of_light: The speed of light in the system.

    Returns:
        Array of time values in the fixed frame of reference.
    """

    t_scale = 1.0 + solver.z_max / (speed_of_light * solver.t_max)
    return solver.tlist * t_scale


@renamed_arg("mb_solve", "solver")
def rabi_freq(
    solver: PropagationSolver,
    field_idx: int,
    speed_of_light: float,
    part: str = "real",
    interp_kind: str = "linear",
) -> np.ndarray:
    """Return the field results shifted to the fixed (lab) frame of reference
    given a speed-of-light by interpolation. Can return the real part or
    abs value.

    Args:
        solver: A PropagationSolver object
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
        rabi_freq_zt = np.abs(solver._Omegas_zt[field_idx])
    elif part == "real":
        rabi_freq_zt = np.real(solver._Omegas_zt[field_idx])
    else:
        raise ValueError('Invalid part. Try "abs" or "real"')

    rabi_freq_intp = RegularGridInterpolator(
        (solver.zlist, solver.tlist),
        rabi_freq_zt,
        method=interp_kind,
        bounds_error=False,
        fill_value=0.0,
    )

    rabi_freq_fixed = np.zeros(solver._Omegas_zt[field_idx].shape, dtype=float)

    for i, z_i in enumerate(solver.zlist):
        t_shifted = t_list(solver, speed_of_light) - z_i / speed_of_light
        pts = np.column_stack([np.full(len(t_shifted), z_i), t_shifted])
        rabi_freq_fixed[i] = rabi_freq_intp(pts)

    return rabi_freq_fixed
