""" Unit tests for the spectral analysis module."""

import os
import unittest

import numpy as np

from maxwellbloch import t_funcs

class TestGaussian(unittest.TestCase):

    def test_areas_pi(self):
        """Test Gaussian areas as multiples of pi.
        """
        FWHM = 0.1
        tlist = np.linspace(0., 1., 201)
        t_func = t_funcs.gaussian(1)
        for n in np.linspace(0.0, 10.0, 11):
            ampl = n*np.sqrt(4.*np.pi*np.log(2)/FWHM**2)/(2*np.pi)  # nπ area
            t_args = {'ampl_1': ampl*2*np.pi, 'fwhm_1': FWHM, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)
            self.assertAlmostEqual(area, n*np.pi, places=3)

    def test_n_pi_areas_pi(self):
        """Test Gaussian areas as multiples of pi given n_pi arg.
        """
        FWHM = 0.1
        tlist = np.linspace(0., 1., 201)
        t_func = t_funcs.gaussian(1)
        for n_pi in np.linspace(0.0, 10.0, 11):
            t_args = {'n_pi_1': n_pi, 'fwhm_1': FWHM, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)
            self.assertAlmostEqual(area, n_pi*np.pi, places=3)

        # TODO: Check error if ampl and n_pi passed
        # TODO: Check error if neither passed.

class TestSech(unittest.TestCase):
    
    def test_areas_pi(self):
        """Test Gaussian areas as multiples of pi.
        """
        SECH_FWHM_CONV = 1./2.6339157938
        FWHM = 0.1
        width = FWHM*SECH_FWHM_CONV # [τ]
        for n in np.linspace(0.0, 10.0, 11):
            ampl = n/width/(2*np.pi) # nπ area
            tlist = np.linspace(0., 1., 201)
            t_func = t_funcs.sech(1)
            t_args = {'ampl_1': ampl*2*np.pi, 'width_1': width, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)
            self.assertAlmostEqual(area, n*np.pi, places=3)
