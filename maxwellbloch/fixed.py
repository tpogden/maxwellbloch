# -*- coding: utf-8 -*-

""" TODO:

"""

import numpy as np
from scipy import interpolate

def t_list(mb_solve, speed_of_light):
    """ Return the time points shifted to the fixed (lab) frame of
        reference given a speed-of-light.

        Args:
            mb_solve: An MBSolve object
            speed_of_light: The speed of light in the system.

        Returns:
            Array of time values in the fixed frame of reference.
    """

    t_scale = 1.0 + mb_solve.z_max/(speed_of_light * mb_solve.t_max)
    return mb_solve.tlist*t_scale

def rabi_freq_abs(mb_solve, field_idx, speed_of_light, interp_kind='linear'):
    """ Return the solved field results shifted to the fixed (lab) frame of
        reference given a speed-of-light by interpolation.

        Args:
            mb_solve: An MBSolve object
            field_idx: The field to return
            speed_of_light: The speed of light in the system
            interp_kind: The kind of spline interpolation to use ('linear',
                'cubic' or 'quintic')

        Returns:
            Array[num_fields, num_space_points, num_time_points] of field
            values in the fixed frame of reference.
    """


    rabi_freq_abs_intp = interpolate.interp2d(mb_solve.tlist, mb_solve.zlist,
        np.abs(mb_solve.Omegas_zt[field_idx]), bounds_error=False,
        fill_value=0., kind=interp_kind)

    rabi_freq_abs_fixed = np.zeros(mb_solve.Omegas_zt[field_idx].shape,
        dtype=np.float)

    for i, z_i in enumerate(mb_solve.zlist):
        rabi_freq_abs_fixed[i] = rabi_freq_abs_intp(
            t_list(mb_solve, speed_of_light) - z_i / speed_of_light, z_i)

    return rabi_freq_abs_fixed

