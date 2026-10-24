"""
Unit tests for the angmom module.

Thomas Ogden <t@ogden.eu>
"""

import unittest

import numpy as np

from clerq import angmom


class TestCalcClebschHF(unittest.TestCase):
    def test_Rb_87_5s12_5p12(self):
        """Testing Clebsch-Gordan coefficients for the Rb87 D1 line hyperfine
            transitions. These coefficients are from Tables 15-20 in ref [1].

        Notes:
            [1]: Daniel A. Steck, 'Rubidium 87 D Line Data' available online
            at http://steck.us/alkalidata. (revision 2.1.4, 23 Dec 2010).
        """

        J_a = 1 / 2
        J_b = 1 / 2
        I_a = 3 / 2
        I_b = 3 / 2

        # Table 15
        def chf_15(F_b, mF_a):
            return angmom.calc_clebsch_hf(
                J_a=J_a,
                I_a=I_a,
                F_a=2,
                mF_a=mF_a,
                J_b=J_b,
                I_b=I_b,
                F_b=F_b,
                mF_b=mF_a + 1,
                q=-1,
            )

        self.assertAlmostEqual(chf_15(F_b=1, mF_a=-2), np.sqrt(1 / 2))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=-1), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=0), np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=1), np.sqrt(0))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=2), np.sqrt(0))

        self.assertAlmostEqual(chf_15(F_b=2, mF_a=-2), np.sqrt(1 / 6))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=-1), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=0), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=1), np.sqrt(1 / 6))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=2), np.sqrt(0))

        # Table 16
        def chf_16(F_b, mF_a):
            return angmom.calc_clebsch_hf(
                J_a=J_a,
                I_a=I_a,
                F_a=2,
                mF_a=mF_a,
                J_b=J_b,
                I_b=I_b,
                F_b=F_b,
                mF_b=mF_a,
                q=0,
            )

        self.assertAlmostEqual(chf_16(F_b=1, mF_a=-2), np.sqrt(0))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=-1), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=0), np.sqrt(1 / 3))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=1), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=2), np.sqrt(0))

        self.assertAlmostEqual(chf_16(F_b=2, mF_a=-2), -np.sqrt(1 / 3))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=-1), -np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=0), np.sqrt(0))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=1), np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=2), np.sqrt(1 / 3))

        # Table 17
        def chf_17(F_b, mF_a):
            return angmom.calc_clebsch_hf(
                J_a=J_a,
                I_a=I_a,
                F_a=2,
                mF_a=mF_a,
                J_b=J_b,
                I_b=I_b,
                F_b=F_b,
                mF_b=mF_a - 1,
                q=1,
            )

        self.assertAlmostEqual(chf_17(F_b=1, mF_a=-2), np.sqrt(0))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=-1), np.sqrt(0))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=0), np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=1), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=2), np.sqrt(1 / 2))

        self.assertAlmostEqual(chf_17(F_b=2, mF_a=-2), np.sqrt(0))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=-1), -np.sqrt(1 / 6))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=0), -np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=1), -np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=2), -np.sqrt(1 / 6))

        # Table 18
        def chf_18(F_b, mF_a):
            return angmom.calc_clebsch_hf(
                J_a=J_a,
                I_a=I_a,
                F_a=1,
                mF_a=mF_a,
                J_b=J_b,
                I_b=I_b,
                F_b=F_b,
                mF_b=mF_a + 1,
                q=-1,
            )

        self.assertAlmostEqual(chf_18(F_b=1, mF_a=-1), -np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_18(F_b=1, mF_a=0), -np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_18(F_b=1, mF_a=1), np.sqrt(0))

        self.assertAlmostEqual(chf_18(F_b=2, mF_a=-1), -np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_18(F_b=2, mF_a=0), -np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_18(F_b=2, mF_a=1), -np.sqrt(1 / 2))

        # Table 19
        def chf_19(F_b, mF_a):
            return angmom.calc_clebsch_hf(
                J_a=J_a,
                I_a=I_a,
                F_a=1,
                mF_a=mF_a,
                J_b=J_b,
                I_b=I_b,
                F_b=F_b,
                mF_b=mF_a,
                q=0,
            )

        self.assertAlmostEqual(chf_19(F_b=1, mF_a=-1), np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_19(F_b=1, mF_a=0), np.sqrt(0))
        self.assertAlmostEqual(chf_19(F_b=1, mF_a=1), -np.sqrt(1 / 12))

        self.assertAlmostEqual(chf_19(F_b=2, mF_a=-1), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_19(F_b=2, mF_a=0), np.sqrt(1 / 3))
        self.assertAlmostEqual(chf_19(F_b=2, mF_a=1), np.sqrt(1 / 4))

        # Table 20
        def chf_20(F_b, mF_a):
            return angmom.calc_clebsch_hf(
                J_a=J_a,
                I_a=I_a,
                F_a=1,
                mF_a=mF_a,
                J_b=J_b,
                I_b=I_b,
                F_b=F_b,
                mF_b=mF_a - 1,
                q=1,
            )

        self.assertAlmostEqual(chf_20(F_b=1, mF_a=-1), np.sqrt(0))
        self.assertAlmostEqual(chf_20(F_b=1, mF_a=0), np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_20(F_b=1, mF_a=1), np.sqrt(1 / 12))

        self.assertAlmostEqual(chf_20(F_b=2, mF_a=-1), -np.sqrt(1 / 2))
        self.assertAlmostEqual(chf_20(F_b=2, mF_a=0), -np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_20(F_b=2, mF_a=1), -np.sqrt(1 / 12))


class TestWigner3j(unittest.TestCase):
    def test_known_values(self):
        """Closed-form values: (j j 0; m -m 0) = (-1)^(j-m) / sqrt(2j+1), etc."""
        w = angmom.wigner_3j
        self.assertAlmostEqual(w(0, 0, 0, 0, 0, 0), 1.0)
        self.assertAlmostEqual(w(1 / 2, 1 / 2, 0, 1 / 2, -1 / 2, 0), np.sqrt(1 / 2))
        self.assertAlmostEqual(w(1 / 2, 1 / 2, 0, -1 / 2, 1 / 2, 0), -np.sqrt(1 / 2))
        self.assertAlmostEqual(w(1, 1, 0, 0, 0, 0), -np.sqrt(1 / 3))
        self.assertAlmostEqual(w(1, 1, 0, 1, -1, 0), np.sqrt(1 / 3))
        self.assertAlmostEqual(w(1, 1, 1, 1, -1, 0), np.sqrt(1 / 6))
        self.assertAlmostEqual(w(1, 1, 2, 0, 0, 0), np.sqrt(2 / 15))
        self.assertAlmostEqual(w(1, 1, 2, 1, -1, 0), np.sqrt(1 / 30))
        self.assertAlmostEqual(w(2, 2, 2, 0, 0, 0), -np.sqrt(2 / 35))

    def test_zero_m_sum_nonzero(self):
        self.assertEqual(angmom.wigner_3j(1, 1, 1, 1, 1, 0), 0)
        self.assertEqual(angmom.wigner_3j(1, 1, 2, 0, 0, 1), 0)

    def test_zero_triangle_violation(self):
        self.assertEqual(angmom.wigner_3j(1, 1, 3, 0, 0, 0), 0)
        self.assertEqual(angmom.wigner_3j(3, 1, 1, 0, 0, 0), 0)
        self.assertEqual(angmom.wigner_3j(1 / 2, 1 / 2, 2, 1 / 2, -1 / 2, 0), 0)

    def test_zero_m_out_of_bounds(self):
        self.assertEqual(angmom.wigner_3j(1, 1, 1, 2, -1, -1), 0)
        self.assertEqual(angmom.wigner_3j(1, 1, 1, -1, 2, -1), 0)
        self.assertEqual(angmom.wigner_3j(1, 1, 1, -1, -1, 2), 0)

    def test_zero_odd_sum_all_m_zero(self):
        """(j1 j2 j3; 0 0 0) vanishes when j1 + j2 + j3 is odd."""
        self.assertAlmostEqual(angmom.wigner_3j(1, 1, 1, 0, 0, 0), 0)
        self.assertAlmostEqual(angmom.wigner_3j(1, 2, 2, 0, 0, 0), 0)
        self.assertAlmostEqual(angmom.wigner_3j(2, 2, 1, 0, 0, 0), 0)

    def test_zero_parity_mismatch(self):
        """j and m must differ by an integer."""
        self.assertEqual(angmom.wigner_3j(1, 1, 1, 1 / 2, -1 / 2, 0), 0)

    def test_non_half_integer_raises(self):
        with self.assertRaises(ValueError):
            angmom.wigner_3j(0.3, 1, 1, 0, 0, 0)

    def test_even_permutation_symmetry(self):
        w = angmom.wigner_3j
        args = (3 / 2, 1, 1 / 2, 1 / 2, -1, 1 / 2)
        j1, j2, j3, m1, m2, m3 = args
        ref = w(*args)
        self.assertAlmostEqual(w(j2, j3, j1, m2, m3, m1), ref)
        self.assertAlmostEqual(w(j3, j1, j2, m3, m1, m2), ref)

    def test_odd_permutation_phase(self):
        w = angmom.wigner_3j
        j1, j2, j3, m1, m2, m3 = (3 / 2, 1, 1 / 2, 1 / 2, -1, 1 / 2)
        phase = (-1) ** (j1 + j2 + j3)
        ref = w(j1, j2, j3, m1, m2, m3)
        self.assertAlmostEqual(w(j2, j1, j3, m2, m1, m3), phase * ref)
        self.assertAlmostEqual(w(j1, j3, j2, m1, m3, m2), phase * ref)

    def test_m_sign_flip_phase(self):
        w = angmom.wigner_3j
        j1, j2, j3, m1, m2, m3 = (3 / 2, 1, 1 / 2, 1 / 2, -1, 1 / 2)
        phase = (-1) ** (j1 + j2 + j3)
        self.assertAlmostEqual(
            w(j1, j2, j3, -m1, -m2, -m3), phase * w(j1, j2, j3, m1, m2, m3)
        )

    def test_orthogonality(self):
        """sum_{m1,m2} (2 j3 + 1) (j1 j2 j3; m1 m2 m3)^2 = 1 for allowed m3."""
        w = angmom.wigner_3j
        for j1, j2, j3 in [(1, 1, 2), (3 / 2, 1, 1 / 2), (2, 3 / 2, 1 / 2)]:
            for m3 in np.arange(-j3, j3 + 1):
                total = sum(
                    (2 * j3 + 1) * w(j1, j2, j3, m1, m3 * -1 - m1, m3) ** 2
                    for m1 in np.arange(-j1, j1 + 1)
                )
                self.assertAlmostEqual(total, 1.0)


class TestWigner6j(unittest.TestCase):
    def test_known_values(self):
        w = angmom.wigner_6j
        self.assertAlmostEqual(w(0, 0, 0, 0, 0, 0), 1.0)
        self.assertAlmostEqual(w(1, 1, 1, 1, 1, 1), 1 / 6)
        self.assertAlmostEqual(w(1 / 2, 1 / 2, 1, 1 / 2, 1 / 2, 1), 1 / 6)
        self.assertAlmostEqual(w(2, 2, 2, 2, 2, 2), -3 / 70)
        # {a b c; 0 c b} = (-1)^(a+b+c) / sqrt((2b+1)(2c+1))
        self.assertAlmostEqual(w(1, 1 / 2, 1 / 2, 0, 1 / 2, 1 / 2), 1 / 2)
        self.assertAlmostEqual(w(1, 1, 1, 0, 1, 1), -1 / 3)
        self.assertAlmostEqual(w(2, 1, 1, 0, 1, 1), 1 / 3)

    def test_zero_triangle_violation(self):
        w = angmom.wigner_6j
        self.assertEqual(w(1, 1, 3, 1, 1, 1), 0)
        self.assertEqual(w(1, 1, 1, 1, 1, 3), 0)
        self.assertEqual(w(1, 1, 1, 3, 1, 1), 0)
        self.assertEqual(w(1, 1, 1, 1, 3, 1), 0)

    def test_zero_non_integer_triad_sum(self):
        self.assertEqual(angmom.wigner_6j(1 / 2, 1, 1, 1, 1, 1), 0)

    def test_non_half_integer_raises(self):
        with self.assertRaises(ValueError):
            angmom.wigner_6j(0.3, 1, 1, 1, 1, 1)

    def test_column_permutation_symmetry(self):
        w = angmom.wigner_6j
        ref = w(3 / 2, 1, 1 / 2, 1, 3 / 2, 1)
        self.assertAlmostEqual(w(1, 3 / 2, 1 / 2, 3 / 2, 1, 1), ref)
        self.assertAlmostEqual(w(1 / 2, 1, 3 / 2, 1, 3 / 2, 1), ref)

    def test_upper_lower_swap_symmetry(self):
        """Swapping upper and lower arguments in any two columns is a symmetry."""
        w = angmom.wigner_6j
        ref = w(3 / 2, 1, 1 / 2, 1, 3 / 2, 1)
        self.assertAlmostEqual(w(1, 3 / 2, 1 / 2, 3 / 2, 1, 1), ref)
        self.assertAlmostEqual(w(1, 1, 1, 3 / 2, 3 / 2, 1 / 2), ref)
