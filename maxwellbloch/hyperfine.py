""" Calculate factors for hyperfine structure in single-electron atoms.
"""

import sys
import numpy as np

class LevelNL(object):
    """ Represents a nL level of a single-electron atom.
    
    Examples:
        Rb_87_5s = LevelNL(n=5, I=1.5, L=1, S=0.5)

    Notes:
        - The J_energies are set _relative_ to the nL level energy.
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
        return "<LevelNL :: %s>" % self.__dict__

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
        return np.arange(abs(self.L - self.S), self.L + self.S + 1)

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
        return "<LevelJ :: %s>" % self.__dict__

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

        return np.arange(abs(self.J - self.I), self.J + self.I + 1)


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
        return "<LevelF :: %s>" % self.__dict__

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

        return np.arange(-self.F, self.F + 1)


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


def main():

    return 0

if __name__ == '__main__':
    STATUS = main()
    sys.exit(STATUS)
