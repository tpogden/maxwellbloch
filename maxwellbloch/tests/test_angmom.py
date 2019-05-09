""" 
Unit tests for the angmom module. 

Thomas Ogden <t@ogden.eu>
"""

import unittest

import numpy as np

from maxwellbloch import angmom

class TestCalcClebschHF(unittest.TestCase):

    def test_Rb_87_5s12_5p12(self):
        """ Testing Clebsch-Gordan coefficients for the Rb87 D1 line hyperfine
            transitions. These coefficients are from Tables 15-20 in ref [1].

        Notes:
            [1]: Daniel A. Steck, 'Rubidium 87 D Line Data' available online
            at http://steck.us/alkalidata. (revision 2.1.4, 23 Dec 2010).
        """

        J_a = 1/2
        J_b = 1/2
        I_a = 3/2
        I_b = 3/2

        # Table 15
        chf_15 = lambda F_b, mF_a: angmom.calc_clebsch_hf(J_a=J_a, I_a=I_a,
            F_a=2, mF_a=mF_a, J_b=J_b, I_b=I_b, F_b=F_b, mF_b=mF_a + 1, q=-1)

        self.assertAlmostEqual(chf_15(F_b=1, mF_a=-2), np.sqrt(1/2))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=-1), np.sqrt(1/4))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=0), np.sqrt(1/12))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=1), np.sqrt(0))
        self.assertAlmostEqual(chf_15(F_b=1, mF_a=2), np.sqrt(0))

        self.assertAlmostEqual(chf_15(F_b=2, mF_a=-2), np.sqrt(1 / 6))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=-1), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=0), np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=1), np.sqrt(1/6))
        self.assertAlmostEqual(chf_15(F_b=2, mF_a=2), np.sqrt(0))

        # Table 16
        chf_16 = lambda F_b, mF_a: angmom.calc_clebsch_hf(J_a=J_a, I_a=I_a,
            F_a=2, mF_a=mF_a, J_b=J_b, I_b=I_b, F_b=F_b, mF_b=mF_a, q=0)

        self.assertAlmostEqual(chf_16(F_b=1, mF_a=-2), np.sqrt(0))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=-1), np.sqrt(1/4))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=0), np.sqrt(1/3))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=1), np.sqrt(1/4))
        self.assertAlmostEqual(chf_16(F_b=1, mF_a=2), np.sqrt(0))

        self.assertAlmostEqual(chf_16(F_b=2, mF_a=-2), -np.sqrt(1/3))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=-1), -np.sqrt(1/12))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=0), np.sqrt(0))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=1), np.sqrt(1/12))
        self.assertAlmostEqual(chf_16(F_b=2, mF_a=2), np.sqrt(1/3))

        # Table 17
        def chf_17(F_b, mF_a): return angmom.calc_clebsch_hf(J_a=J_a, I_a=I_a,
            F_a=2, mF_a=mF_a, J_b=J_b, I_b=I_b, F_b=F_b, mF_b=mF_a - 1, q=1)

        self.assertAlmostEqual(chf_17(F_b=1, mF_a=-2), np.sqrt(0))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=-1), np.sqrt(0))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=0), np.sqrt(1/12))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=1), np.sqrt(1/4))
        self.assertAlmostEqual(chf_17(F_b=1, mF_a=2), np.sqrt(1/2))

        self.assertAlmostEqual(chf_17(F_b=2, mF_a=-2), np.sqrt(0))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=-1), -np.sqrt(1/6))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=0), -np.sqrt(1/4))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=1), -np.sqrt(1/4))
        self.assertAlmostEqual(chf_17(F_b=2, mF_a=2), -np.sqrt(1/6))

        # Table 18
        chf_18 = lambda F_b, mF_a: angmom.calc_clebsch_hf(J_a=J_a, I_a=I_a,
            F_a=1, mF_a=mF_a, J_b=J_b, I_b=I_b, F_b=F_b, mF_b=mF_a+1, q=-1)

        self.assertAlmostEqual(chf_18(F_b=1, mF_a=-1), -np.sqrt(1/12))
        self.assertAlmostEqual(chf_18(F_b=1, mF_a=0), -np.sqrt(1/12))
        self.assertAlmostEqual(chf_18(F_b=1, mF_a=1), np.sqrt(0))

        self.assertAlmostEqual(chf_18(F_b=2, mF_a=-1), -np.sqrt(1/12))
        self.assertAlmostEqual(chf_18(F_b=2, mF_a=0), -np.sqrt(1/4))
        self.assertAlmostEqual(chf_18(F_b=2, mF_a=1), -np.sqrt(1/2))

        # Table 19
        chf_19 = lambda F_b, mF_a: angmom.calc_clebsch_hf(J_a=J_a, I_a=I_a,
            F_a=1, mF_a=mF_a, J_b=J_b, I_b=I_b, F_b=F_b, mF_b=mF_a, q=0)

        self.assertAlmostEqual(chf_19(F_b=1, mF_a=-1), np.sqrt(1/12))
        self.assertAlmostEqual(chf_19(F_b=1, mF_a=0), np.sqrt(0))
        self.assertAlmostEqual(chf_19(F_b=1, mF_a=1), -np.sqrt(1/12))

        self.assertAlmostEqual(chf_19(F_b=2, mF_a=-1), np.sqrt(1/4))
        self.assertAlmostEqual(chf_19(F_b=2, mF_a=0), np.sqrt(1/3))
        self.assertAlmostEqual(chf_19(F_b=2, mF_a=1), np.sqrt(1/4))

        # Table 20
        def chf_20(F_b, mF_a): return angmom.calc_clebsch_hf(J_a=J_a, I_a=I_a,
            F_a=1, mF_a=mF_a, J_b=J_b, I_b=I_b, F_b=F_b, mF_b=mF_a - 1, q=1)

        self.assertAlmostEqual(chf_20(F_b=1, mF_a=-1), np.sqrt(0))
        self.assertAlmostEqual(chf_20(F_b=1, mF_a=0), np.sqrt(1 / 12))
        self.assertAlmostEqual(chf_20(F_b=1, mF_a=1), np.sqrt(1 / 12))

        self.assertAlmostEqual(chf_20(F_b=2, mF_a=-1), -np.sqrt(1 / 2))
        self.assertAlmostEqual(chf_20(F_b=2, mF_a=0), -np.sqrt(1 / 4))
        self.assertAlmostEqual(chf_20(F_b=2, mF_a=1), -np.sqrt(1 / 12))
