""" Helper functions for calculating macroscopic properties of atomic 
    ensembles.
"""

from scipy import constants as c

a0 = c.physical_constants['atomic unit of length'][0]
e = c.physical_constants['atomic unit of charge'][0]
hbar = c.physical_constants['Planck constant over 2 pi'][0]
eps0 = c.physical_constants['electric constant'][0]

def convert_dipole_moment_atomic_to_si(dipole_moment_ea0):
    """ Convert a dipole moment from atomic units [e*a_0] to SI
        units [C m].
    """
    return dipole_moment_ea0*e*a0

def propagation_coefficient_si(dipole_moment, wavenumber):
    """ Propagation coefficient, g, in SI units
        Args:
            dipole_moment [C m]
            wavenumber [m^-1]
        Returns:
            propagation_coefficient (g) [m^2 s^-1]
        References:
            Thesis Eqn (2.53)
    """
    return dipole_moment**2*wavenumber/(2*eps0*hbar)

def propagation_coefficient_MHz_cm2(dipole_moment, wavenumber):
    """
        Args:
            dipole_moment [C m]
            wavenumber [m^-1]
        Returns:
            propagation_coefficient (g) [2π MHz cm^2]
    """
    Hz_m2_to_2piMHz_cm2 = 1e-6*1e4/(2*c.pi)
    
    return propagation_coefficient_si(dipole_moment, 
        wavenumber)*Hz_m2_to_2piMHz_cm2

def interaction_strength_si(dipole_moment, wavenumber, number_density):
    """ Interaction strength, N*g, in SI units
        Args:
            dipole_moment [C m]
            wavenumber [m^-1]
            number_density [m^-3]
        Returns:
            interaction_strength (N*g) [s^-1 m^-1]
    """            
    return propagation_coefficient_si(dipole_moment, wavenumber)*number_density

def interaction_strength_MHz_per_cm(dipole_moment, wavenumber, number_density):
    """    
        Args:
            dipole_moment [C m]
            wavenumber [m^-1]
            number_density [cm^-3]
         Returns:
            interaction_strength (N*g) [2π MHz cm^-1]
    """

    g_MHz_cm2 = propagation_coefficient_MHz_cm2(dipole_moment, wavenumber)    
    return number_density*g_MHz_cm2 # [2π MHz cm^-1]
