# -*- coding: utf-8 -*-

"""Conversions between laboratory (experimental) units and clerq's units.

clerq works in a consistent set of "natural" units: frequencies are cyclic
(``rabi_freq``, ``decays`` rates, detunings and ``thermal_width`` are ``Omega /
2 pi``, ``Gamma / 2 pi`` and so on), and time is the inverse. The helpers here
fix that unit system to

- frequencies in MHz and times in microseconds,
- lengths in cm and number densities in cm^-3,
- wavelengths in nm,
- intensities in mW/cm^2, powers in mW and beam sizes in cm.

They return values ready to pass to a :class:`~clerq.System` or
:class:`~clerq.PropagationSolver`.

Two-level conventions: the transition has natural (spontaneous) linewidth
``decay_rate`` = ``Gamma / 2 pi`` and wavelength ``wavelength``, and
``factor`` scales the dipole moment squared relative to the closed cycling
transition (e.g. a squared Clebsch-Gordan coefficient). With ``factor=1`` the
transition is a cycling one with ``Gamma = omega^3 d^2 / (3 pi eps0 hbar
c^3)``.
"""

from __future__ import annotations

import numpy as np
from scipy import constants as _const

_NM = 1e-9
_CM = 1e-2
_MHZ = 1e6
# mW/cm^2 -> W/m^2
_MW_CM2 = 1e-3 / 1e-4

AMU_KG = _const.atomic_mass
"""Atomic mass unit, in kg."""


def saturation_intensity(decay_rate: float, wavelength: float) -> float:
    """Saturation intensity of a cycling transition, ``pi h c Gamma / (3 lambda^3)``.

    Args:
        decay_rate: Natural linewidth ``Gamma / 2 pi``, in MHz.
        wavelength: Transition wavelength, in nm.

    Returns:
        ``I_sat`` in mW/cm^2. For the Rb87 D2 line (6.065 MHz, 780.24 nm) this
        is 1.67 mW/cm^2.
    """
    gamma = 2 * np.pi * decay_rate * _MHZ
    i_sat = np.pi * _const.h * _const.c * gamma / (3 * (wavelength * _NM) ** 3)
    return i_sat / _MW_CM2


def intensity_to_rabi_freq(
    intensity: float, decay_rate: float, wavelength: float, factor: float = 1.0
) -> float:
    """Cyclic Rabi frequency of a field of the given intensity.

    ``Omega = Gamma * sqrt(factor * I / (2 I_sat))``, returned as ``Omega / 2
    pi``.

    Args:
        intensity: Intensity at the atoms, in mW/cm^2.
        decay_rate: Natural linewidth ``Gamma / 2 pi``, in MHz.
        wavelength: Transition wavelength, in nm.
        factor: Dipole moment squared relative to the cycling transition.

    Returns:
        ``Omega / 2 pi`` in MHz, the value for a clerq ``rabi_freq``.
    """
    i_sat = saturation_intensity(decay_rate, wavelength)
    return decay_rate * np.sqrt(factor * np.asarray(intensity) / (2 * i_sat))


def rabi_freq_to_intensity(
    rabi_freq: float, decay_rate: float, wavelength: float, factor: float = 1.0
) -> float:
    """Intensity, in mW/cm^2, for a cyclic Rabi frequency ``Omega / 2 pi`` in MHz.

    The inverse of :func:`intensity_to_rabi_freq`.
    """
    i_sat = saturation_intensity(decay_rate, wavelength)
    return 2 * i_sat * (np.asarray(rabi_freq) / decay_rate) ** 2 / factor


def gaussian_beam_peak_intensity(power: float, waist: float) -> float:
    """Peak intensity of a circular Gaussian beam, ``2 P / (pi w^2)``.

    Args:
        power: Total beam power, in mW.
        waist: The ``1/e^2`` intensity radius, in cm.

    Returns:
        Peak intensity in mW/cm^2.
    """
    return 2 * power / (np.pi * waist**2)


def interaction_strength(
    decay_rate: float, wavelength: float, factor: float = 1.0
) -> float:
    """Interaction strength per atom, ``g / 2 pi`` with ``g = 3 lambda^2 Gamma / (8 pi)``.

    The propagation equation is ``dOmega/dz = 2 i N g rho``, so multiplying
    by the number density ``N`` (cm^-3) gives the clerq
    ``interaction_strengths`` entry, in MHz/cm, for a density profile of
    amplitude 1 (the default). A weak resonant field is then absorbed at
    ``N * factor * 3 lambda^2 / (2 pi)`` per cm in intensity.

    Args:
        decay_rate: Natural linewidth ``Gamma / 2 pi``, in MHz.
        wavelength: Transition wavelength, in nm.
        factor: Dipole moment squared relative to the cycling transition.

    Returns:
        ``g / 2 pi`` in MHz cm^2.
    """
    lam_cm = wavelength * _NM / _CM
    return 3 * lam_cm**2 * decay_rate * factor / (8 * np.pi)


def resonant_optical_depth(
    number_density: float, length: float, wavelength: float, factor: float = 1.0
) -> float:
    """Resonant intensity optical depth ``N sigma_0 L`` of a cold, weak-field medium.

    Args:
        number_density: Atoms per cm^3.
        length: Medium length, in cm.
        wavelength: Transition wavelength, in nm.
        factor: Dipole moment squared relative to the cycling transition.

    Returns:
        The optical depth ``OD``; intensity transmission is ``exp(-OD)``. Doppler
        broadening reduces the peak.
    """
    lam_cm = wavelength * _NM / _CM
    sigma0 = factor * 3 * lam_cm**2 / (2 * np.pi)
    return number_density * sigma0 * length


def thermal_width(temperature: float, mass: float, wavelength: float) -> float:
    """The ``thermal_width`` of ``velocity_classes`` for a vapor.

    The Doppler weights in clerq are ``exp(-(delta / w)^2)`` with ``w = u /
    lambda`` and ``u = sqrt(2 k_B T / m)`` the most probable speed.

    Args:
        temperature: Temperature, in kelvin.
        mass: Atomic mass, in amu.
        wavelength: Transition wavelength, in nm.

    Returns:
        ``w`` in MHz.
    """
    u = np.sqrt(2 * _const.k * temperature / (mass * AMU_KG))
    return u / (wavelength * _NM) / _MHZ
