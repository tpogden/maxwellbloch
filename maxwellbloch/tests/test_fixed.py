# -*- coding: utf-8 -*-

""" Unit tests for the fixed module for returning MBSolve results in the fixed
frame of reference.
"""

import os

from maxwellbloch import mb_solve, fixed

import numpy as np

import unittest

# Absolute path of tests/json directory, so that tests can be called from
# different directories.
JSON_DIR = os.path.abspath(os.path.join(__file__, '../', 'json'))

class TestTlist(unittest.TestCase):
    """ Unit tests of the tlist method. """

    def test_tlist_max(self):
        """ When the speed of light is 1.5, for t_max = 1.0 and z_max = 1.0,
            the t_max_fixed should be 1 + 1/1.5 = 5/3. """

        mb_solve_00 = mb_solve.MBSolve()

        speed_of_light = 1.5  # speed of light
        t_max_fixed_expected = mb_solve_00.t_max + \
            1 / (speed_of_light * mb_solve_00.z_max)

        t_max_fixed = fixed.t_list(mb_solve_00, speed_of_light)[-1]

        self.assertAlmostEqual(t_max_fixed, t_max_fixed_expected)

class TestRabiFreqAbs(unittest.TestCase):
    """ Unite tests of the rabi_freq_abs method. """

    def test_rabi_freq_abs(self):

        json_path = os.path.join(JSON_DIR, "mb_solve_no_atoms.json")
        mb_solve_01 = mb_solve.MBSolve().from_json(json_path)

        mb_solve_01.mbsolve(recalc=False)

        speed_of_light = 0.5

        t_list_fixed = fixed.t_list(mb_solve_01, speed_of_light)
        rabi_freq_abs_fixed = fixed.rabi_freq_abs(mb_solve_01, 0,
            speed_of_light, interp_kind='linear')

        # The max of the the field in both frames should be the same
        self.assertAlmostEqual(np.max(rabi_freq_abs_fixed),
            np.max(mb_solve_01.Omegas_zt), places=0)
