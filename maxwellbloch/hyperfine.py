""" Calculate factors for hyperfine structure in single-electron atoms.
"""

import sys
import numpy as np


class LevelJ(object):
    """ Represents a J fine structure level. """

    def __init__(self, J, energy):
        self.J = J
        self.F_levels = []
        self.energy = energy

    def __repr__(self):
        return "<LevelJ :: %s>" % self.__dict__

    def build_F_levels(self, I, F_energies, mF_energies=None):
        for i, F in enumerate(self.get_F_range(I)):
            self.add_F_level(LevelF(F, F_energies[i]))
            # TODO: I think this check should be outside loop
            if (mF_energies == None):
                self.F_levels[i].build_mF_levels()
            elif (len(mF_energies) == self.get_num_F_levels(I)):
                self.F_levels[i].build_mF_levels(mF_energies[i])
            else:
                pass  # should raise an error here

    def add_F_level(self, F_level):
        self.F_levels.append(F_level)

    def get_F_range(self, I):
        return np.arange(abs(self.J - I), self.J + I + 1)

    def get_num_F_levels(self, I):
        return self.get_F_range(I).size

class LevelF(object):
    """ Represents an F hyperfine structure level. 
    
        Attributes:
            F (float): Total atomic angular momentum number F.
            energy (float): Energy of the level.
            mf_levels (list, length 2F+1): Energies of 2F+1 hyperfine
                sublevels.

        Notes:
            The magnitude of F can take values in the range 
            `|J - I| <= F <= J + 1`.
    """

    def __init__(self, F, energy, mF_energies=None):
        """ Inits LevelF.

        Args:
            F (float): Total atomic angular momentum number F.
            energy (float): Energy of the level.
            mf_energies (list(LevelMF), length 2F+1): List of 2F+1 hyperfine
                sublevels.
        """
        self.F = F
        self.energy = energy
        self.mF_levels = []
        
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

        mF_range = self.get_mF_range()
        if not mF_energies:
            mF_energies = [self.energy for i in mF_range]
        if len(mF_energies) != len(mF_range):
            raise ValueError("mF_energies is not the correct length.")
        for i, mF in enumerate(mF_range):
            self.mF_levels.append(LevelMF(mF, mF_energies[i]))

        return self.mF_levels

    def get_mF_range(self):
        """ Returns a range representing the m_F angular momentum sublevels
            [-F, -F+1, ..., F-1, F].
        """
        
        return np.arange(-self.F, self.F + 1)

    # def get_num_mF_levels(self):

    #     return self.get_mF_range().size

class LevelMF(object):
    """ Represents an m_F hyperfine sublevel. 
    
    Attributes:
        mF (float): Angular momentum number m_F
        energy (float): The energy of the level
    """

    def __init__(self, mF, energy):
        self.mF = mF
        self.energy = energy

    def __repr__(self):
        return "<LevelMF :: %s>" % self.__dict__


def main():

    my_F = LevelF(F=2.0, energy=0.0)
    my_F.build_mF_levels()
    print(my_F)

    return 0

if __name__ == '__main__':
    STATUS = main()
    sys.exit(STATUS)
