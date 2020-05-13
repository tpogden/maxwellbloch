"""
Unit tests for the utility module.

Thomas Ogden <t@ogden.eu>
"""

import os
import unittest

import numpy as np

from maxwellbloch import utility, t_funcs

class TestHalfMaxRoots(unittest.TestCase):

    def test_gaussian(self):
        """ Tests a Gaussian pulse for the correct half max and width. """
        t_list = tlist = np.linspace(0., 1., 101)
        y_list = t_funcs.gaussian(1)(tlist, 
            args={'ampl_1':1.0, 'fwhm_1':0.1, 'centre_1':0.5})
        half_max, r1, r2 = utility.half_max_roots(t_list, y_list)
        self.assertAlmostEqual(half_max, 0.5)
        self.assertAlmostEqual(r1, 0.45)
        self.assertAlmostEqual(r2, 0.55)

class TestFullWidthAtHalfMax(unittest.TestCase):

    def test_gaussian(self):
        """ Tests a Gaussian pulse for the correct FWHM. """
        t_list = tlist = np.linspace(0., 1., 101)
        y_list = t_funcs.gaussian(1)(tlist, 
            args={'ampl_1':1.0, 'fwhm_1':0.1, 'centre_1':0.5})
        self.assertAlmostEqual(utility.full_width_at_half_max(t_list, y_list),
            0.1)
