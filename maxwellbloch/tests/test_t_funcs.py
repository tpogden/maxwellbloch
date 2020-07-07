""" Unit tests for the spectral analysis module."""

import os
import unittest

import numpy as np

from maxwellbloch import t_funcs, utility

class TestGaussian(unittest.TestCase):

    def test_areas_pi(self):
        """Test Gaussian areas as multiples of pi.
        """
        FWHM = 0.1
        tlist = np.linspace(0., 1., 201)
        t_func = t_funcs.gaussian(1)
        for n in np.linspace(1.0, 10.0, 10):
            ampl = n*np.sqrt(4.*np.pi*np.log(2)/FWHM**2)/(2*np.pi)  # nπ area
            t_args = {'ampl_1': ampl, 'fwhm_1': FWHM, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)*2*np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, 
                t_func(tlist, t_args))
            self.assertAlmostEqual(area, n*np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_areas_pi_n_pi(self):
        """Test Gaussian areas as multiples of pi given n_pi arg.
        """
        FWHM = 0.1
        tlist = np.linspace(0., 1., 201)
        t_func = t_funcs.gaussian(1)
        for n_pi in np.linspace(1.0, 10.0, 10):
            t_args = {'n_pi_1': n_pi, 'fwhm_1': FWHM, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)*2*np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, 
                t_func(tlist, t_args))
            self.assertAlmostEqual(area, n_pi*np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_ampl_and_n_pi(self):
        """Test that KeyError is raised if both ampl and n_pi args set.
        """
        tlist = np.linspace(0., 1., 201)
        t_args = {'n_pi_1': 2.0, 'ampl_1': 1.0, 'fwhm_1': 0.1, 'centre_1': 0.5}
        t_func = t_funcs.gaussian(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)

    def test_no_ampl_nor_n_pi(self):
        tlist = np.linspace(0., 1., 201)
        t_args = {'fwhm_1': 0.1, 'centre_1': 0.5}
        t_func = t_funcs.gaussian(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)

class TestSech(unittest.TestCase):
    
    def test_areas_pi(self):
        """Test sech areas as multiples of pi.
        """
        SECH_FWHM_CONV = 1./2.6339157938
        FWHM = 0.1
        width = FWHM*SECH_FWHM_CONV # [τ]
        tlist = np.linspace(0., 1., 201)
        t_func = t_funcs.sech(1)
        for n in np.linspace(1.0, 10.0, 10):
            ampl = n/width/(2*np.pi) # nπ area
            t_args = {'ampl_1': ampl, 'width_1': width, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)*2*np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, 
                t_func(tlist, t_args))
            self.assertAlmostEqual(area, n*np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_areas_pi_n_pi(self):
        """Test sech areas as multiples of pi given n_pi arg.
        """
        SECH_FWHM_CONV = 1./2.6339157938
        FWHM = 0.1
        width = FWHM*SECH_FWHM_CONV # [τ]
        tlist = np.linspace(0., 1., 201)
        t_func = t_funcs.sech(1)
        for n in np.linspace(1.0, 10.0, 10):
            t_args = {'n_pi_1': n, 'width_1': width, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)*2*np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, 
                t_func(tlist, t_args))
            self.assertAlmostEqual(area, n*np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_areas_pi_n_pi_fwhm(self):
        """Test sech areas as multiples of pi given n_pi and fwhm args.
        """
        tlist = np.linspace(0., 1., 201)
        t_func = t_funcs.sech(1)
        FWHM = 0.1
        for n in np.linspace(1.0, 10.0, 10):
            t_args = {'n_pi_1': n, 'fwhm_1': FWHM, 'centre_1': 0.5}
            area = np.trapz(t_func(tlist, t_args), tlist)*2*np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, 
                t_func(tlist, t_args))
            self.assertAlmostEqual(area, n*np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)
    # TODO: Test the FWHM is correct

    def test_ampl_and_n_pi(self):
        """Test that KeyError is raised if both ampl and n_pi args set.
        """
        tlist = np.linspace(0., 1., 201)
        t_args = {'n_pi_1': 2.0, 'ampl_1': 1.0, 'width_1': 0.1, 'centre_1': 0.5}
        t_func = t_funcs.sech(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)

    def test_no_ampl_nor_n_pi(self):
        tlist = np.linspace(0., 1., 201)
        t_args = {'width_1': 0.1, 'centre_1': 0.5}
        t_func = t_funcs.sech(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)
