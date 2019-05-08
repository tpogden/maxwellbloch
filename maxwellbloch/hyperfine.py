""" Calculate factors for hyperfine structure in single-electron atoms.
"""

import sys
import json
from itertools import product

from maxwellbloch.angmom import calc_clebsch_hf

import numpy as np

class Atom1e(object):
    """ Represents an atom object and contains a list of J levels. """

    def __init__(self, element=None, isotope=None):

        self.element = element
        self.isotope = isotope
        self.J_levels = []

    def __repr__(self):

        return self.to_json_str()

    def add_J_level(self, J_level):

        self.J_levels.append(J_level)

    def get_J_level_idx_map(self):
        """ Maps each element of the mF_list to a J_level index. """

        J_level_idx_map = []
        for i, J_level in enumerate(self.J_levels):
            for F_level in J_level.F_levels:
                for _ in F_level.mF_levels:
                    J_level_idx_map.append(i)
        return J_level_idx_map

    def get_mF_list(self):
        """ Unnests the mF levels into a single list of dicts. """

        mF_list = []

        for J_level in self.J_levels:
            item_I = J_level.I
            item_J = J_level.J
            for F_level in J_level.F_levels:
                item_F = F_level.F
                for mF_level in F_level.mF_levels:
                    item_mF = mF_level.mF
                    item_energy = mF_level.energy
                    mF_dict = {'I':item_I, 'J':item_J, 'F':item_F, 
                        'mF':item_mF, 'energy':item_energy}
                    mF_list.append(mF_dict)

        return mF_list

    def get_num_mF_levels(self):
        """ Returns the total number of mF levels in all J levels. """

        return len(self.get_mF_list())

    def get_energies(self):
        """ Returns a list of energies for each mF level in all J levels. """

        return [mF_level['energy'] for mF_level in self.get_mF_list()]

    def get_coupled_levels(self, J_level_idx_a, J_level_idx_b):
        """ Returns a list of pairs of mF level indexes, representing all
            pairs of mF levels between two J levels. """

        J_level_idx_map = self.get_J_level_idx_map()
        a_levels = [i for i, idx in enumerate(J_level_idx_map) 
            if idx == J_level_idx_a]
        b_levels = [i for i, idx in enumerate(J_level_idx_map)
            if idx == J_level_idx_b]

        # product returns iterator of tuples, convert to list of lists
        return [list(i) for i in product(a_levels, b_levels)]

    def get_clebsch_hf_factors(self, J_level_idx_a, J_level_idx_b, q):
        """ Returns a list of Clebsch-Gordan coefficients for the hyperfine 
            transition dipole matrix elements for each coupled level pair.
        
        Args:
            J_level_idx_a (int): J level the transition is from (lower level)
            J_level_idx_b (int): J level the transition is to (upper level)
            q (int): The field polarisation. Choose from [-1, 0, 1].

        Returns:
            (list): factors, length of mF_list

        """

        factors = []
        mF_list = self.get_mF_list()
        coupled_levels = self.get_coupled_levels(J_level_idx_a, J_level_idx_b)
        for cl in coupled_levels:
            a = mF_list[cl[0]]
            b = mF_list[cl[1]]
            clebsch_hf = calc_clebsch_hf(J_a=a['J'], I_a=a['I'], F_a=a['F'],
                mF_a=a['mF'], J_b=b['J'], I_b=b['I'], F_b=b['F'], mF_b=b['mF'], 
                q=q)
            factors.append(clebsch_hf)
        return factors

    # TODO: iso_factors (mix of all 3 polarisations)

    def get_decay_factors(self, J_level_idx_a, J_level_idx_b):
        """ Returns a list of factors for the collapse operators for each 
            hyperfine coupled level pair.

            Args: 
                J_level_idx_a(int): J level the transition is from (lower level)
                J_level_idx_b(int): J level the transition is to(upper level)
                q(int): The field polarisation. Choose from [-1, 0, 1].

            Returns: (list): factors, length of mF_list
        """

        decay_factors = []
        mF_list = self.get_mF_list()
        coupled_levels = self.get_coupled_levels(J_level_idx_a, J_level_idx_b)
        for cl in coupled_levels:
            a = mF_list[cl[0]]
            b = mF_list[cl[1]]
            clebsch_hf = 0.0
            # Only one of these should be nonzero so OK to sum them
            for q in [-1, 0, 1]:
                clebsch_hf += calc_clebsch_hf(J_a=a['J'], I_a=a['I'], 
                    F_a=a['F'], mF_a=a['mF'], J_b=b['J'], I_b=b['I'], 
                    F_b=b['F'], mF_b=b['mF'], q=q)
            decay_factors.append(clebsch_hf)
        return decay_factors

    def to_json_str(self):
        """ Return a JSON string representation of the LevelJ object.

        # TODO: This could be a decorator as we use it for all classes.

        Returns:
            (string) JSON representation of the LevelJ object.
        """

        return json.dumps(self.get_json_dict(), indent=2, separators=None, 
            sort_keys=True)

    def get_json_dict(self):

        json_dict = {"element": self.element,
                     "isotope": self.isotope,
                     "J_levels": [i.get_json_dict() for i in self.J_levels]}
        return json_dict

class LevelNL(object):
    """ Represents a nL level of a single-electron atom.
    
    Examples:
        Rb_87_5s = LevelNL(n=5, I=1.5, L=1, S=0.5)

    Notes:
        - The J_energies are set _relative_ to the nL level energy.
    TODO: I may not need this anymore, now Atom1e takes a list of J_levels
    """

    def __init__(self, n, I, L, S, energy=0.0, J_energies=None, 
        F_energies=None, mF_energies=None):
        """
        Args:
            n (float): principal atomic number.
            I (float): Nuclear spin atomic number.
            L (float): Orbital angular momenutm number.
            S (float): Spin angular momentum number.
            energy (float): Energy of the nL level
            J_energies (list of float): List of energies of each J level.
            F_energies (list of list of float): List of energies of the F 
                levels, relative to each J level.
            mF_energies (list of list of list of float): List of list of lists 
                containing mF_energies for each F level within each J level. 
                The length of each sublist must be 2F+1.
        """

        self.n = n
        self.I = I 
        self.L = L
        self.S = S
        self.energy = energy
        self.J_levels = self.build_J_levels(J_energies, F_energies, 
            mF_energies)

    def __repr__(self):
        return self.to_json_str() #"<LevelNL :: %s>" % self.__dict__

    def build_J_levels(self, J_energies=None, F_energies=None, 
        mF_energies=None):
        
        self.J_levels = []
        J_range = self.get_J_range()

        if not J_energies:
            J_energies = [0.0 for i in J_range]
        if not F_energies:
            F_energies = [None for i in J_range]
        if not mF_energies:
            mF_energies = [None for i in J_range]

        if len(J_energies) != len(J_range):
            raise ValueError("J_energies is not the correct length.")
        if len(F_energies) != len(J_range):
            raise ValueError("F_energies is not the correct length.")
        if len(mF_energies) != len(J_range):
            raise ValueError("mF_energies is not the correct length.")

        for i, J in enumerate(J_range):
            self.J_levels.append(LevelJ(self.I, J, self.energy + J_energies[i], 
                F_energies[i], mF_energies[i]))

        return self.J_levels

    def get_J_range(self):
        return np.arange(abs(self.L - self.S), self.L + self.S + 1, 
            dtype=float)

    def to_json_str(self):
        """ Return a JSON string representation of the LevelJ object.

        # TODO: This could be a decorator as we use it for all classes.

        Returns:
            (string) JSON representation of the LevelJ object.
        """

        return json.dumps(self.get_json_dict(), indent=2, separators=None, 
            sort_keys=True)

    def get_json_dict(self):

        json_dict = {"n": self.n,
                     "I": self.I,
                     "L": self.L,
                     "S": self.S,
                     "energy": self.energy,
                     "J_levels": [i.get_json_dict() for i in self.J_levels]}
        return json_dict

class LevelJ(object):
    """ Represents a J fine structure level and holds its hyperfine structure
        sublevels. 

    Examples:
        Rb_87_5p12 = LevelJ(I=1.5, J=0.5)

    Notes:
        - The magnitude of J can take values in the range 
            `|L - S| <= J <= L + S`.
        - The F levels take values in the range `|J - I| <= F <= J + I`.
        - The F_energies are set _relative_ to the J level energy.
    """

    def __init__(self, I, J, energy=0.0, F_energies=None, mF_energies=None):
        """
        Args:
            I (float): Nuclear spin atomic number.
            J (float): Orbital angular momentum number.
            energy (float): Energy of the J level
            F_energies (list of float): List of energies of the F levels, 
                relative to the J level.
            mF_energies (list of list of float): List of lists containing 
                mF_energies for each F level. The length of each list must be 
                2F+1.

        Notes:
            - I and J must be integer or half-integer.
        """

        if ((2*I != round(2*I)) | (2*J != round(2*J))):
            raise ValueError('I and J must be integers or half-integers.')

        self.J = J
        self.I = I
        self.energy = energy
        self.build_F_levels(F_energies, mF_energies)

    def __repr__(self):
        return self.to_json_str() #"<LevelJ :: %s>" % self.__dict__

    def build_F_levels(self, F_energies=None, mF_energies=None):
        """ Builds the hyperfine structure F levels of the J level.

        Args:
            F_energies (list of float): List of energies of the F levels.
            mF_energies (list of list of float): List of lists containing 
                mF_energies for each F level. The length of each list must be 
                2F+1.

        """
        self.F_levels = []
        F_range = self.get_F_range()

        if not F_energies:
            F_energies = [0.0 for i in F_range]
        if not mF_energies:
            mF_energies = [None for i in F_range]
        if len(F_energies) != len(F_range):
            raise ValueError("F_energies is not the correct length.")
        if len(mF_energies) != len(F_range):
            raise ValueError("mF_energies is not the correct length.")
        for i, F in enumerate(F_range):
            self.F_levels.append(LevelF(F, self.energy + F_energies[i], 
                mF_energies[i]))

        return self.F_levels

    def get_F_range(self):
        """ The range of the F levels is given by `|J - I| <= F <= J + I`. """ 

        return np.arange(abs(self.J - self.I), self.J + self.I + 1, 
            dtype=float)

    def get_json_dict(self):

        json_dict = {"I": self.I,
                     "J": self.J,
                     "energy": self.energy,
                     "F_levels": [F.get_json_dict() for F in self.F_levels]}
        return json_dict

    def to_json_str(self):
        """ Return a JSON string representation of the LevelJ object.

        # TODO: This could be a decorator as we use it for all classes.

        Returns:
            (string) JSON representation of the LevelJ object.
        """

        return json.dumps(self.get_json_dict(), indent=2, separators=None, 
            sort_keys=True)

class LevelF(object):
    """ Represents an F hyperfine structure level and holds its magnetic 
        sublevels mF. 
    
    Attributes:
        F (float): Total atomic angular momentum number F.
        energy (float): Energy of the level.
        mf_levels (list, length 2F+1): Energies of 2F+1 hyperfine
            sublevels.

    Notes:
        - The magnitude of F can take values in the range 
            `|J - I| <= F <= J + I`.
        - The mF_energies are set _relative_ to the F level energy.
    """

    def __init__(self, F, energy=0.0, mF_energies=None):
        """ 
        Args:
            F (float): Total atomic angular momentum number F.
            energy (float): Energy of the level.
            mf_energies (list(LevelMF), length 2F+1): List of 2F+1 hyperfine
            sublevels.
        """
        
        self.F = F
        self.energy = energy
        self.build_mF_levels(mF_energies)

    def __repr__(self):
        return self.to_json_str() #"<LevelF :: %s>" % self.__dict__

    def build_mF_levels(self, mF_energies=None):
        """ Builds the mF sublevels of the F level. 
        
        Args:
            mF_energies  (list, length 2F+1): Energies of the 2F+1 hyperfine
                sublevels.
        Returns:
            self.mF_levels
        """
        self.mF_levels = []
        mF_range = self.get_mF_range()
        if not mF_energies:
            mF_energies = [0.0 for i in mF_range]
        if len(mF_energies) != len(mF_range):
            raise ValueError("mF_energies is not the correct length.")
        for i, mF in enumerate(mF_range):
            self.mF_levels.append(LevelMF(mF=mF, energy=self.energy + 
                mF_energies[i]))

        return self.mF_levels

    def get_mF_range(self):
        """ Returns a range representing the m_F angular momentum sublevels
            [-F, -F+1, ..., F-1, F].
        """

        return np.arange(-self.F, self.F + 1, dtype=float)

    def get_json_dict(self):

        json_dict = {"F": self.F,
                     "energy": self.energy,
                     "mF_levels": [mF.__dict__ for mF in self.mF_levels]}
        return json_dict

    def to_json_str(self):
        """ Return a JSON string representation of the LevelF object.

        Returns:
            (string) JSON representation of the LevelF object.
        """

        return json.dumps(self.get_json_dict(), indent=2, separators=None, 
            sort_keys=True)


class LevelMF(object):
    """ Represents an m_F hyperfine sublevel. 
    
    Attributes:
        mF (float): Angular momentum number m_F
        energy (float): The energy of the level
    """

    def __init__(self, mF, energy=0.0):
        self.mF = mF
        self.energy = energy

    def __repr__(self):
        return "<LevelMF :: %s>" % self.__dict__

    def to_json_str(self):
        """ Return a JSON string representation of the LevelMF object.

        Returns:
            (string) JSON representation of the LevelMF object.
        """


        return json.dumps(self.__dict__, indent=2, separators=None, 
            sort_keys=True)


def main():

    # Rb_87_5s = LevelNL(n=5, I=1.5, L=0, S=0.5)

    # Rb_87_5p = LevelNL(n=5, I=1.5, L=1, S=0.5)

    # print(Rb_87_5s)

    Rb87_5s12 = LevelJ(I=1.5, J=0.5)

    Rb87_5p12 = LevelJ(I=1.5, J=0.5)
    Rb87_5p32 = LevelJ(I=1.5, J=1.5)

    Rb87_5s_5p = Atom1e(element='Rb', isotope='87')
    # atom.add_nL_level(Rb_87_5s)
    # atom.add_nL_level(Rb_87_5p)

    Rb87_5s_5p.add_J_level(Rb87_5s12)
    Rb87_5s_5p.add_J_level(Rb87_5p12)
    # Rb87_5s_5p.add_J_level(Rb87_5p32)

    print(Rb87_5s_5p)

    print(Rb87_5s_5p.get_mF_list())

    print(Rb87_5s_5p.get_J_level_idx_map())

    print(Rb87_5s_5p.get_coupled_levels(0,1))

    print(Rb87_5s_5p.get_energies())

    print('Transition factors:')

    print(Rb87_5s_5p.get_clebsch_hf_factors(0, 1, q=-0))

    print('Decay factors:')

    print(Rb87_5s_5p.get_decay_factors(0, 1))

    return 0

if __name__ == '__main__':
    STATUS = main()
    sys.exit(STATUS)
