""" 
Unit tests for the hyperfine module. 

Thomas Ogden <t@ogden.eu>
"""

import unittest

import numpy as np

from maxwellbloch import hyperfine

        
class TestAtom1eAddFLevel(unittest.TestCase):
    """ Unit tests of the Atom1e.add_F_level method. """

    def test_Rb87_5s12_5p12(self):
        """ TODO: put a levels diagram here. """

        Rb87_5s12_5p12 = hyperfine.Atom1e(element='Rb', isotope='87')
        # Create F levels and add to atom
        Rb87_5s12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5s12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)
        Rb87_5p12_F1 = hyperfine.LevelF(I=1.5, J=1.5, F=1)
        Rb87_5p12_F2 = hyperfine.LevelF(I=1.5, J=1.5, F=2)
        Rb87_5s12_5p12.add_F_level(Rb87_5s12_F1)
        Rb87_5s12_5p12.add_F_level(Rb87_5s12_F2)
        Rb87_5s12_5p12.add_F_level(Rb87_5p12_F1)
        Rb87_5s12_5p12.add_F_level(Rb87_5p12_F2)

        self.assertEqual(Rb87_5s12_5p12.get_num_mF_levels(), 16)

        map = [0]*3 + [1]*5 + [2]*3 + [3]*5
        self.assertEqual(Rb87_5s12_5p12.get_F_level_idx_map(), map)

        self.assertEqual(len(Rb87_5s12_5p12.get_coupled_levels([0], [2])), 3*3)
        self.assertEqual(len(Rb87_5s12_5p12.get_coupled_levels([0], [3])), 3*5)

        self.assertEqual(len(Rb87_5s12_5p12.get_coupled_levels([0], [2, 3])), 
            3*8)
        self.assertEqual(len(Rb87_5s12_5p12.get_coupled_levels([0, 1], [2])), 
            8*3)

class TestAtom1eGetClebschHFFactors(unittest.TestCase):

    def setup_method(self, method):

        Rb87_5s12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5s12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)
        Rb87_5p12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5p12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)

        self.Rb87_5s12_5p12 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F2)

        Rb87_5p32_F0 = hyperfine.LevelF(I=1.5, J=1.5, F=0)
        Rb87_5p32_F1 = hyperfine.LevelF(I=1.5, J=1.5, F=1)
        Rb87_5p32_F2 = hyperfine.LevelF(I=1.5, J=1.5, F=2)        
        Rb87_5p32_F3 = hyperfine.LevelF(I=1.5, J=1.5, F=3)

        self.Rb87_5s12_5p32 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F0)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F3)

    def test_Rb87_5s12_5p12(self):
        """
        References:
            [0]: https://steck.us/alkalidata/rubidium87numbers.pdf
        """
        # self.assertEqual(self.Rb87_5s12_5p12.get_num_mF_levels(), 16)

        Rb87_5s12_F_level_idxs = (0, 1)
        Rb87_5p12_F_level_idxs = (2, 3)

        cl = self.Rb87_5s12_5p12.get_coupled_levels(Rb87_5s12_F_level_idxs, 
            Rb87_5p12_F_level_idxs)

        facts_qm1 = self.Rb87_5s12_5p12.get_clebsch_hf_factors(
            Rb87_5s12_F_level_idxs, Rb87_5p12_F_level_idxs, q=-1)
        facts_q0 = self.Rb87_5s12_5p12.get_clebsch_hf_factors(
            Rb87_5s12_F_level_idxs, Rb87_5p12_F_level_idxs, q=0)
        facts_qp1 = self.Rb87_5s12_5p12.get_clebsch_hf_factors(
            Rb87_5s12_F_level_idxs, Rb87_5p12_F_level_idxs, q=1)

        # For each upper mF level, get the sum of all couplings mF' 
        # Eqn (40) in Ref [0]. These should sum to (2J + 1)/(2J' + 1) = 1
        # "The interpretation of this symmetry is simply that all the excited 
        # state sublevels decay at the same rate \Gamma, and the decaying 
        # population “branches” into various ground state sublevels."
        for upper_mF_level in range(8, 16):
            coupled = [upper_mF_level in i for i in cl]

            factor_sq_sum = (np.sum(facts_qm1[coupled]**2) + 
                np.sum(facts_q0[coupled]**2) + np.sum(facts_qp1[coupled]**2))

            self.assertAlmostEqual(factor_sq_sum, 1.0) # (2J+1)/(2J'+1) = 1

class TestAtom1eGetStrengthFactor(unittest.TestCase):

    def setup_method(self, method):

        Rb87_5s12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5s12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)
        Rb87_5p12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5p12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)

        self.Rb87_5s12_5p12 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F2)

        Rb87_5p32_F0 = hyperfine.LevelF(I=1.5, J=1.5, F=0)
        Rb87_5p32_F1 = hyperfine.LevelF(I=1.5, J=1.5, F=1)
        Rb87_5p32_F2 = hyperfine.LevelF(I=1.5, J=1.5, F=2)        
        Rb87_5p32_F3 = hyperfine.LevelF(I=1.5, J=1.5, F=3)

        self.Rb87_5s12_5p32 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F0)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F3)

    def test_Rb87_5s12_5p12_strength_factors(self):
        """ Test strength factors against known values from Steck [0] Table 8.
    
        Refs:
            [0]: https://steck.us/alkalidata/rubidium87numbers.pdf
        """

        # The loops are to check each lower state, should be the same for each
        for i in range(2):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p12.get_strength_factor(0, 2, i), 1.0/6) # F1 -> F1
        for i in range(2):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p12.get_strength_factor(0, 3, i), 5.0/6) # F1 -> F2
        for i in range(4):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p12.get_strength_factor(1, 2, i), 1.0/2) # F2 -> F1
        for i in range(4):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p12.get_strength_factor(1, 3, i), 1.0/2) # F2 -> F2

        # The sum of S_{FF'} over upper F levels should be 1.
        for j in range(2):
            self.assertAlmostEqual(
                self.Rb87_5s12_5p12.get_strength_factor(j, 2) + 
                self.Rb87_5s12_5p12.get_strength_factor(j, 3), 1.0)

    def test_Rb87_5s12_5p32_strength_factors(self):
        """ Test strength factors against known values from Steck [0] Table 8.
    
        Refs:
            [0]: https://steck.us/alkalidata/rubidium87numbers.pdf
        """

        # The loops are to check each lower state, should be the same for each
        for i in range(2):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(0, 2, i), 1./6) # F1 -> F0
        for i in range(2):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(0, 3, i), 5./12) # F1 -> F1
        for i in range(2):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(0, 4, i), 5./12) # F1 -> F2
        for i in range(2):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(0, 5, i), 0.0) # F1 -> F3

        for i in range(4):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(1, 2, i), 0.0) # F1 -> F0
        for i in range(4):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(1, 3, i), 1./20) # F1 -> F1
        for i in range(4):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(1, 4, i), 1./4) # F1 -> F2
        for i in range(4):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p32.get_strength_factor(1, 5, i), 7./10) # F1 -> F3

        # The sum of S_{FF'} over upper F levels should be 1.
        for j in range(2):
            self.assertAlmostEqual(
            self.Rb87_5s12_5p12.get_strength_factor(j, 2) + 
            self.Rb87_5s12_5p12.get_strength_factor(j, 3) +
            self.Rb87_5s12_5p12.get_strength_factor(j, 4) +
            self.Rb87_5s12_5p12.get_strength_factor(j, 5), 1.0)

    def test_Rb87_5s12_5p32(self):

        Rb87_5s12_F_level_idxs = (0, 1)
        Rb87_5p32_F_level_idxs = (2, 3, 4, 5)
        cl = self.Rb87_5s12_5p32.get_coupled_levels(Rb87_5s12_F_level_idxs, 
            Rb87_5p32_F_level_idxs)

        facts_qm1 = self.Rb87_5s12_5p32.get_clebsch_hf_factors(
            Rb87_5s12_F_level_idxs, Rb87_5p32_F_level_idxs, q=-1)
        facts_q0 = self.Rb87_5s12_5p32.get_clebsch_hf_factors(
            Rb87_5s12_F_level_idxs, Rb87_5p32_F_level_idxs, q=0)
        facts_qp1 = self.Rb87_5s12_5p32.get_clebsch_hf_factors(
            Rb87_5s12_F_level_idxs, Rb87_5p32_F_level_idxs, q=1)

        # For each upper mF level, get the sum of all couplings mF' 
        # Eqn (40) in Ref [0]. These should sum to (2J + 1)/(2J' + 1) = 0.5
        # "The interpretation of this symmetry is simply that all the excited 
        # state sublevels decay at the same rate \Gamma, and the decaying 
        # population “branches” into various ground state sublevels."
        for upper_mF_level in range(8, 24):
            coupled = [upper_mF_level in i for i in cl]

            factor_sq_sum = (np.sum(facts_qm1[coupled]**2) + 
                np.sum(facts_q0[coupled]**2) + np.sum(facts_qp1[coupled]**2))

            self.assertAlmostEqual(factor_sq_sum, 0.5) # (2J+1)/(2J'+1) = 0.5

class TestAtom1eGetDecayFactors(unittest.TestCase):

    def setup_method(self, method):

        Rb87_5s12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5s12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)
        Rb87_5p12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5p12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)
        
        self.Rb87_5s12_5p12 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F2)

        Rb87_5p32_F0 = hyperfine.LevelF(I=1.5, J=1.5, F=0)
        Rb87_5p32_F1 = hyperfine.LevelF(I=1.5, J=1.5, F=1)
        Rb87_5p32_F2 = hyperfine.LevelF(I=1.5, J=1.5, F=2)        
        Rb87_5p32_F3 = hyperfine.LevelF(I=1.5, J=1.5, F=3)

        self.Rb87_5s12_5p32 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F0)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F3)

    def test_Rb87_5s12_5p12(self):

        Rb87_5s12_F_level_idxs = (0, 1)
        Rb87_5p12_F_level_idxs = (2, 3)
        cl = self.Rb87_5s12_5p12.get_coupled_levels(Rb87_5s12_F_level_idxs, 
            Rb87_5p12_F_level_idxs)
        df = self.Rb87_5s12_5p12.get_decay_factors(Rb87_5s12_F_level_idxs, 
            Rb87_5p12_F_level_idxs)

        # For each upper mF level, get the sum of all couplings mF' 
        # Eqn (40) in Ref [0]. These should sum to (2J + 1)/(2J' + 1) = 1
        # "The interpretation of this symmetry is simply that all the excited 
        # state sublevels decay at the same rate \Gamma, and the decaying 
        # population “branches” into various ground state sublevels."
        for upper_mF_level in range(8, 16):
            coupled = [upper_mF_level in i for i in cl]
            factor_sq_sum = np.sum(df[coupled]**2)
            self.assertAlmostEqual(factor_sq_sum, 1.0) # (2J+1)/(2J'+1) = 1

    def test_Rb87_5s12_5p32(self):

        Rb87_5s12_F_level_idxs = (0, 1)
        Rb87_5p32_F_level_idxs = (2, 3, 4, 5)
        cl = self.Rb87_5s12_5p32.get_coupled_levels(Rb87_5s12_F_level_idxs, 
            Rb87_5p32_F_level_idxs)
        df = self.Rb87_5s12_5p32.get_decay_factors(Rb87_5s12_F_level_idxs, 
            Rb87_5p32_F_level_idxs)

        for upper_mF_level in range(8, 24):
            coupled = [upper_mF_level in i for i in cl]
            factor_sq_sum = np.sum(df[coupled]**2)
            self.assertAlmostEqual(factor_sq_sum, 0.5) # (2J+1)/(2J'+1) = 1/2

class GetClebschHFFactorsIso(unittest.TestCase):

    def setup_method(self, method):

        Rb87_5s12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5s12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)
        Rb87_5p12_F1 = hyperfine.LevelF(I=1.5, J=0.5, F=1)
        Rb87_5p12_F2 = hyperfine.LevelF(I=1.5, J=0.5, F=2)
        
        self.Rb87_5s12_5p12 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F1)
        self.Rb87_5s12_5p12.add_F_level(Rb87_5p12_F2)

        Rb87_5p32_F0 = hyperfine.LevelF(I=1.5, J=1.5, F=0)
        Rb87_5p32_F1 = hyperfine.LevelF(I=1.5, J=1.5, F=1)
        Rb87_5p32_F2 = hyperfine.LevelF(I=1.5, J=1.5, F=2)        
        Rb87_5p32_F3 = hyperfine.LevelF(I=1.5, J=1.5, F=3)

        self.Rb87_5s12_5p32 = hyperfine.Atom1e(element='Rb', isotope='87')
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5s12_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F0)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F1)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F2)
        self.Rb87_5s12_5p32.add_F_level(Rb87_5p32_F3)

    def test_Rb87_5s12_5p12(self):

        Rb87_5s12_F_level_idxs = (0, 1)
        Rb87_5p12_F_level_idxs = (2, 3)

        cl = self.Rb87_5s12_5p12.get_coupled_levels(
            Rb87_5s12_F_level_idxs, Rb87_5p12_F_level_idxs)
        cf = self.Rb87_5s12_5p12.get_clebsch_hf_factors_iso(
            Rb87_5s12_F_level_idxs, Rb87_5p12_F_level_idxs)

        for upper_mF_level in range(8, 16):
            coupled = [upper_mF_level in i for i in cl]
            factor_sq_sum = np.sum(cf[coupled]**2)
            self.assertAlmostEqual(factor_sq_sum, 1.0/3)

    def test_Rb87_5s12_5p32(self):

        Rb87_5s12_F_level_idxs = (0, 1)
        Rb87_5p32_F_level_idxs = (2, 3, 4, 5)

        cl = self.Rb87_5s12_5p32.get_coupled_levels(
            Rb87_5s12_F_level_idxs, Rb87_5p32_F_level_idxs)
        cf = self.Rb87_5s12_5p32.get_clebsch_hf_factors_iso(
            Rb87_5s12_F_level_idxs, Rb87_5p32_F_level_idxs)

        for upper_mF_level in range(8, 24):
            coupled = [upper_mF_level in i for i in cl]
            factor_sq_sum = np.sum(cf[coupled]**2)
            self.assertAlmostEqual(factor_sq_sum, 0.5/3) # (2J+1)/(2J'+1) = 1/2

# class TestLevelNLInit(unittest.TestCase):
#     """ Unit tests of the LevelNL.__init__ method. """

#     def test_Rb_87_5s(self):
#         """ Test the structure built for the Rubidium 87 5s level. """

#         Rb_87_5s = hyperfine.LevelNL(n=5, I=1.5, L=1, S=0.5)
#         self.assertEqual(len(Rb_87_5s.J_levels), 2)

#         # Fine structure
#         self.assertEqual(len(Rb_87_5s.J_levels[0].F_levels), 2)
#         self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels), 4)

#         # Hyperfine structure
#         self.assertEqual(len(Rb_87_5s.J_levels[0].F_levels[0].mF_levels), 3)
#         self.assertEqual(len(Rb_87_5s.J_levels[0].F_levels[1].mF_levels), 5)

#         self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[0].mF_levels), 1)
#         self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[1].mF_levels), 3)
#         self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[2].mF_levels), 5)
#         self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[3].mF_levels), 7)

# class TestLevelJInit(unittest.TestCase):
#     """ Unit tests of the LevelJ.__init__ method. """

#     def test_not_half_ints(self):
#         """ Test that a ValueError is raised if I and J are not both half 
#             integer. """

#         with self.assertRaises(ValueError):
#             hyperfine.LevelJ(I=1.6, J=0.5, energy=0.0)
#         with self.assertRaises(ValueError):
#             hyperfine.LevelJ(I=1.5, J=0.6, energy=0.0)
         
#     def test_no_F_energies(self):
#         """ Test with no setting of F_energies or mF_energies (so should all be 
#             zero) """

#         I = 1.5 # Rb87 nuclear spin
#         J = 0.5 # 5p_{1/2}
#         Rb_87_5p12 = hyperfine.LevelJ(I=I, J=J, energy=0.0, 
#             F_energies=None, mF_energies=None)
#         self.assertEqual(len(Rb_87_5p12.F_levels), 2)
#         self.assertEqual(len(Rb_87_5p12.F_levels[0].mF_levels), 2*1+1)
#         self.assertEqual(len(Rb_87_5p12.F_levels[1].mF_levels), 2*2+1)

#     def test_Rb_87_5p12(self):
#         """ Test building the 5p_{1/2} level of Rubidium 87. """

#         I = 1.5 # Rb87 nuclear spin
#         J = 0.5 # 5p_{1/2}
#         Rb_87_5p12 = hyperfine.LevelJ(I=I, J=J, energy=0.0, 
#             F_energies=[10.0, 20.0])

#         # F numbers are in the range |J - I| <= F <= J + I, so in this case
#         # F=1 and F=2
#         self.assertEqual(len(Rb_87_5p12.F_levels), 2)

#         # Each F level should have 2*F+1 sublevels.
#         self.assertEqual(len(Rb_87_5p12.F_levels[0].mF_levels), 2*1+1)
#         self.assertEqual(len(Rb_87_5p12.F_levels[1].mF_levels), 2*2+1)

class TestLevelFInit(unittest.TestCase):
    """ Unit tests of the LevelF.__init__ method. """

    def test_f_2(self):
        """ Tests that without specified energies, the mF sublevel energies
            are set to zero. """
        
        F = 2
        lf = hyperfine.LevelF(I=1.5, J=0.5, F=F, energy=0.0)
        self.assertEqual(len(lf.mF_levels), 2*F+1)
        for i in lf.mF_levels:
            self.assertEqual(i.energy, 0.0)

    def test_f_2_energies(self):
        """ Tests that the mF sublevel energies are set. """

        F = 2
        mF_energies = [10.0, 20.0, 30.0, 40.0, 50.0]
        lf = hyperfine.LevelF(I=1.5, J=0.5, F=F, energy=0.0, 
            mF_energies=mF_energies)

        self.assertEqual(len(lf.mF_levels), 2*F+1)
        for i, e in enumerate(mF_energies):
            self.assertEqual(e, lf.mF_levels[i].energy)
