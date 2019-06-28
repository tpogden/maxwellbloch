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

    def test_Rb87_5s12_5p32(self):

        Rb87_5s12_F_level_idxs = (0, 1)
        Rb87_5p32_F_level_idxs = (2, 3, 4, 5)

        cl = self.Rb87_5s12_5p32.get_coupled_levels(Rb87_5s12_F_level_idxs, 
            Rb87_5p32_F_level_idxs)

        print(cl)

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
