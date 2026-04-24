# -*- coding: utf-8 -*-

"""Spectral analysis of MBSolve results."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from maxwellbloch.mb_solve import MBSolve


def freq_list(mb_solve: MBSolve) -> np.ndarray:
    """Fourier transform of the tlist into the frequency domain for
        spectral analysis.

    Args:
        mb_solve: An MBSolve object.

        Returns:
            Array[num_time_points] of frequency values.

    """

    t_step = mb_solve.tlist[1] - mb_solve.tlist[0]
    f_list = np.fft.fftfreq(len(mb_solve.tlist), t_step)  # FFT Freq
    return np.fft.fftshift(f_list)


def rabi_freq(mb_solve: MBSolve, field_idx: int) -> np.ndarray:
    """Fourier transform of the field result of field index.

    Args:
        mb_solve: An MBSolve object.
        field_idx: Field to return FFT of.

    Returns:
        Array[num_z_steps, num_t_steps] Field result in frequency domain.

    """

    rabi_freq_zt = mb_solve.Omegas_zt[field_idx]

    return np.fft.fftshift(np.fft.fft(rabi_freq_zt, axis=1), axes=1)


def absorption(mb_solve: MBSolve, field_idx: int, z_idx: int = -1) -> np.ndarray:
    """Field absorption in the frequency domain.

    Args:
        mb_solve: An MBSolve object.
        field_idx: Field to return spectrum of.
        z_idx: z step at which to return absorption.

    Returns:
        Array[num_freq_points] of absorption values.

    Note:
        In the linear regime this is the imaginary part of the linear
        susceptibility (with a factor k/2).
        See TP Ogden thesis Eqn (2.58)
    """

    rabi_freq_abs = np.abs(rabi_freq(mb_solve, field_idx))

    return -np.log(rabi_freq_abs[z_idx] / rabi_freq_abs[0])


def dispersion(mb_solve: MBSolve, field_idx: int, z_idx: int = -1) -> np.ndarray:
    """Field dispersion in the frequency domain.

    Args:
        mb_solve: An MBSolve object.
        field_idx: Field to return spectrum of.
        z_idx: z step at which to return absorption.

    Note:
        In the linear regime this is the real part of the linear
        susceptibility.

        See TP Ogden Thesis Eqn (2.59)
    """

    Omega_freq_angle = np.angle(rabi_freq(mb_solve, field_idx))

    return Omega_freq_angle[0] - Omega_freq_angle[z_idx]


def susceptibility_two_linear_known(
    freq_list: np.ndarray, interaction_strength: float, decay_rate: float
) -> np.ndarray:
    """In the linear regime for a two-level system, the suecpetibility is
    known analytically. This is here for useful comparison, as good
    agreement between a simulated weak field in a two-level system tells us
    that the model is accurate in the linear regime, which gives us
    confidence in the scheme for going beyond the weak field limit.

    Notes:
        See TP Ogden Thesis Eqn(2.61)
    """

    return 1j * interaction_strength / (decay_rate / 2 - 1j * freq_list)


def absorption_two_linear_known(
    freq_list: np.ndarray, interaction_strength: float, decay_rate: float
) -> np.ndarray:
    """The absorption is half the imaginary part of the susecptibility."""

    return (
        susceptibility_two_linear_known(
            freq_list, interaction_strength, decay_rate
        ).imag
        / 2.0
    )


def dispersion_two_linear_known(
    freq_list: np.ndarray, interaction_strength: float, decay_rate: float
) -> np.ndarray:
    """The dispersion is half the real part of the susecptibility."""

    return (
        susceptibility_two_linear_known(
            freq_list, interaction_strength, decay_rate
        ).real
        / 2.0
    )


def voigt_two_linear_known(
    freq_list: np.ndarray, decay_rate: float, thermal_width: float
) -> np.ndarray:
    """Returns the Voigt profile for a two-level system in the linear regime.

    The Voigt profile is the convolution of a Lorentzian with a Gaussian, and
    describes the absorption lineshape for a thermal two-level system.

    Args:
        freq_list: List of frequency detunings from resonance (in 2pi Gamma).
        decay_rate: Spontaneous decay rate of the transition (in 2pi Gamma).
        thermal_width: Width of the lineshape in the same units as decay rate
            (in 2pi Gamma).

    Notes:
        See my thesis section 2.5.6 for more information.
    """
    from scipy.special import wofz

    a = decay_rate / (2 * np.pi * thermal_width)
    b = freq_list / (2 * np.pi * thermal_width)
    s = 1.0j * (0.5 * np.sqrt(np.pi) / (2 * np.pi * thermal_width)) * wofz(b + 0.5j * a)
    return s
