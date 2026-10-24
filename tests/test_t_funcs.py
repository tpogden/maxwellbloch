"""Unit tests for the spectral analysis module."""

import os
import unittest

import numpy as np

from clerq import t_funcs, utility


class TestGaussian(unittest.TestCase):
    def test_areas_pi(self):
        """Test Gaussian areas as multiples of pi."""
        FWHM = 0.1
        tlist = np.linspace(0.0, 1.0, 201)
        t_func = t_funcs.gaussian(1)
        for n in np.linspace(1.0, 10.0, 10):
            ampl = (
                n * np.sqrt(4.0 * np.pi * np.log(2) / FWHM**2) / (2 * np.pi)
            )  # nπ area
            t_args = {"ampl_1": ampl, "fwhm_1": FWHM, "centre_1": 0.5}
            area = np.trapezoid(t_func(tlist, t_args), tlist) * 2 * np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, t_func(tlist, t_args))
            self.assertAlmostEqual(area, n * np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_areas_pi_n_pi(self):
        """Test Gaussian areas as multiples of pi given n_pi arg."""
        FWHM = 0.1
        tlist = np.linspace(0.0, 1.0, 201)
        t_func = t_funcs.gaussian(1)
        for n_pi in np.linspace(1.0, 10.0, 10):
            t_args = {"n_pi_1": n_pi, "fwhm_1": FWHM, "centre_1": 0.5}
            area = np.trapezoid(t_func(tlist, t_args), tlist) * 2 * np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, t_func(tlist, t_args))
            self.assertAlmostEqual(area, n_pi * np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_ampl_and_n_pi(self):
        """Test that KeyError is raised if both ampl and n_pi args set."""
        tlist = np.linspace(0.0, 1.0, 201)
        t_args = {"n_pi_1": 2.0, "ampl_1": 1.0, "fwhm_1": 0.1, "centre_1": 0.5}
        t_func = t_funcs.gaussian(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)

    def test_no_ampl_nor_n_pi(self):
        tlist = np.linspace(0.0, 1.0, 201)
        t_args = {"fwhm_1": 0.1, "centre_1": 0.5}
        t_func = t_funcs.gaussian(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)


class TestSech(unittest.TestCase):
    def test_areas_pi(self):
        """Test sech areas as multiples of pi."""
        SECH_FWHM_CONV = 1.0 / 2.6339157938
        FWHM = 0.1
        width = FWHM * SECH_FWHM_CONV  # [τ]
        tlist = np.linspace(0.0, 1.0, 201)
        t_func = t_funcs.sech(1)
        for n in np.linspace(1.0, 10.0, 10):
            ampl = n / width / (2 * np.pi)  # nπ area
            t_args = {"ampl_1": ampl, "width_1": width, "centre_1": 0.5}
            area = np.trapezoid(t_func(tlist, t_args), tlist) * 2 * np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, t_func(tlist, t_args))
            self.assertAlmostEqual(area, n * np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_areas_pi_n_pi(self):
        """Test sech areas as multiples of pi given n_pi arg."""
        SECH_FWHM_CONV = 1.0 / 2.6339157938
        FWHM = 0.1
        width = FWHM * SECH_FWHM_CONV  # [τ]
        tlist = np.linspace(0.0, 1.0, 201)
        t_func = t_funcs.sech(1)
        for n in np.linspace(1.0, 10.0, 10):
            t_args = {"n_pi_1": n, "width_1": width, "centre_1": 0.5}
            area = np.trapezoid(t_func(tlist, t_args), tlist) * 2 * np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, t_func(tlist, t_args))
            self.assertAlmostEqual(area, n * np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_areas_pi_n_pi_fwhm(self):
        """Test sech areas as multiples of pi given n_pi and fwhm args."""
        tlist = np.linspace(0.0, 1.0, 201)
        t_func = t_funcs.sech(1)
        FWHM = 0.1
        for n in np.linspace(1.0, 10.0, 10):
            t_args = {"n_pi_1": n, "fwhm_1": FWHM, "centre_1": 0.5}
            area = np.trapezoid(t_func(tlist, t_args), tlist) * 2 * np.pi
            fwhm_test = utility.full_width_at_half_max(tlist, t_func(tlist, t_args))
            self.assertAlmostEqual(area, n * np.pi, places=3)
            self.assertAlmostEqual(fwhm_test, FWHM)

    def test_fwhm(self):
        """sech with fwhm arg should have the stated FWHM."""
        tlist = np.linspace(0.0, 1.0, 401)
        t_func = t_funcs.sech(1)
        for fwhm in [0.05, 0.1, 0.2]:
            t_args = {"n_pi_1": 1.0, "fwhm_1": fwhm, "centre_1": 0.5}
            fwhm_measured = utility.full_width_at_half_max(tlist, t_func(tlist, t_args))
            self.assertAlmostEqual(fwhm_measured, fwhm, places=3)

    def test_ampl_and_n_pi(self):
        """Test that KeyError is raised if both ampl and n_pi args set."""
        tlist = np.linspace(0.0, 1.0, 201)
        t_args = {"n_pi_1": 2.0, "ampl_1": 1.0, "width_1": 0.1, "centre_1": 0.5}
        t_func = t_funcs.sech(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)

    def test_no_ampl_nor_n_pi(self):
        tlist = np.linspace(0.0, 1.0, 201)
        t_args = {"width_1": 0.1, "centre_1": 0.5}
        t_func = t_funcs.sech(1)
        with self.assertRaises(KeyError):
            t_func(tlist, t_args)


class TestSin2(unittest.TestCase):
    def test_shape(self):
        """sin2 peaks at the midpoint, vanishes at and outside the edges."""
        tlist = np.linspace(0.0, 1.0, 1001)
        t_args = {"ampl_1": 2.0, "on_1": 0.2, "off_1": 0.8}
        y = t_funcs.sin2(1)(tlist, t_args)
        self.assertAlmostEqual(y.max(), 2.0)
        self.assertAlmostEqual(tlist[np.argmax(y)], 0.5)
        self.assertTrue(np.all(y[tlist < 0.2] == 0.0))
        self.assertTrue(np.all(y[tlist > 0.8] == 0.0))

    def test_fwhm_and_area(self):
        tlist = np.linspace(0.0, 1.0, 2001)
        t_args = {"ampl_1": 1.0, "on_1": 0.2, "off_1": 0.8}
        y = t_funcs.sin2(1)(tlist, t_args)
        self.assertAlmostEqual(utility.full_width_at_half_max(tlist, y), 0.3, places=3)
        self.assertAlmostEqual(np.trapezoid(y, tlist), 0.3, places=4)


class TestSuperGaussian(unittest.TestCase):
    def test_fwhm_any_order(self):
        """The fwhm arg is the full width at half maximum for every order."""
        tlist = np.linspace(0.0, 1.0, 2001)
        for order in [1, 2, 5, 10]:
            t_args = {"ampl_1": 1.0, "fwhm_1": 0.2, "centre_1": 0.5, "order_1": order}
            y = t_funcs.super_gaussian(1)(tlist, t_args)
            # The spline-based utility fails on the sharp edges of high orders.
            fwhm = np.sum(y >= 0.5) * (tlist[1] - tlist[0])
            self.assertAlmostEqual(fwhm, 0.2, places=2)

    def test_order_1_is_gaussian(self):
        tlist = np.linspace(0.0, 1.0, 201)
        sg = t_funcs.super_gaussian(1)(
            tlist, {"ampl_1": 1.0, "fwhm_1": 0.1, "centre_1": 0.5, "order_1": 1}
        )
        g = t_funcs.gaussian(1)(tlist, {"ampl_1": 1.0, "fwhm_1": 0.1, "centre_1": 0.5})
        np.testing.assert_allclose(sg, g, atol=1e-12)

    def test_flat_top(self):
        """A high order is flat at the centre and has vanished outside the width."""
        t_args = {"ampl_1": 1.0, "fwhm_1": 0.4, "centre_1": 0.5, "order_1": 10}
        f = t_funcs.super_gaussian(1)
        self.assertAlmostEqual(f(0.5 + 0.1, t_args), 1.0, places=4)
        self.assertAlmostEqual(f(0.5 + 0.3, t_args), 0.0, places=4)


class TestExpDecayRise(unittest.TestCase):
    def test_decay(self):
        tlist = np.linspace(0.0, 1.0, 201)
        t_args = {"ampl_1": 1.5, "on_1": 0.3, "tau_1": 0.1}
        y = t_funcs.exp_decay(1)(tlist, t_args)
        self.assertTrue(np.all(y[tlist < 0.3] == 0.0))
        self.assertAlmostEqual(t_funcs.exp_decay(1)(0.3, t_args), 1.5)
        self.assertAlmostEqual(t_funcs.exp_decay(1)(0.4, t_args), 1.5 / np.e)

    def test_rise_is_time_reverse(self):
        tlist = np.linspace(0.0, 1.0, 201)
        t_args = {"ampl_1": 1.5, "off_1": 0.7, "tau_1": 0.1}
        y = t_funcs.exp_rise(1)(tlist, t_args)
        self.assertTrue(np.all(y[tlist > 0.7] == 0.0))
        self.assertAlmostEqual(t_funcs.exp_rise(1)(0.7, t_args), 1.5)
        self.assertAlmostEqual(t_funcs.exp_rise(1)(0.6, t_args), 1.5 / np.e)

    def test_no_overflow(self):
        """The masked-out side must not produce inf or nan."""
        t_args = {"ampl_1": 1.0, "off_1": 0.0, "tau_1": 1e-3}
        y = t_funcs.exp_rise(1)(np.linspace(0.0, 10.0, 11), t_args)
        self.assertTrue(np.all(np.isfinite(y)))


class TestSinc(unittest.TestCase):
    def test_centre(self):
        """sinc peaks at the centre arg, and at 0 when it is absent."""
        tlist = np.linspace(-1.0, 1.0, 2001)
        f = t_funcs.sinc(1)
        y0 = f(tlist, {"ampl_1": 1.0, "width_1": 5.0})
        yc = f(tlist, {"ampl_1": 1.0, "width_1": 5.0, "centre_1": 0.3})
        self.assertAlmostEqual(tlist[np.argmax(y0)], 0.0)
        self.assertAlmostEqual(tlist[np.argmax(yc)], 0.3, places=3)
