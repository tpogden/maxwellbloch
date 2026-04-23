# -*- coding: utf-8 -*-

"""Unit tests for OBBase via its OBAtom subclass.

OBBase is not intended to be instantiated directly (num_states=0 makes its
methods meaningless), so we exercise its methods through a minimal OBAtom
instance, which is the direct subclass of OBBase.
"""

import unittest

import numpy as np
import qutip as qu

from maxwellbloch import ob_atom


def make_two_state():
    """Return a minimal two-state OBAtom."""
    return ob_atom.OBAtom(num_states=2)


class TestSigma(unittest.TestCase):
    """Tests for OBBase.sigma()."""

    def test_sigma_ground_to_excited(self):
        """sigma(0, 1) should be |0><1|, i.e. the raising operator."""
        obs = make_two_state()
        sig = obs.sigma(0, 1)
        self.assertIsInstance(sig, qu.Qobj)
        expected = np.array([[0, 1], [0, 0]], dtype=complex)
        np.testing.assert_array_equal(sig.full(), expected)

    def test_sigma_excited_to_ground(self):
        """sigma(1, 0) should be |1><0|, i.e. the lowering operator."""
        obs = make_two_state()
        sig = obs.sigma(1, 0)
        expected = np.array([[0, 0], [1, 0]], dtype=complex)
        np.testing.assert_array_equal(sig.full(), expected)

    def test_sigma_diagonal(self):
        """sigma(0, 0) should be |0><0|, the ground-state projector."""
        obs = make_two_state()
        sig = obs.sigma(0, 0)
        expected = np.array([[1, 0], [0, 0]], dtype=complex)
        np.testing.assert_array_equal(sig.full(), expected)


class TestSetH0(unittest.TestCase):
    """Tests for OBBase.set_H_0()."""

    def test_empty_energies_gives_zero_hamiltonian(self):
        """Empty energy list → all-zero bare Hamiltonian."""
        obs = make_two_state()
        H_0 = obs.set_H_0([])
        self.assertIsInstance(H_0, qu.Qobj)
        np.testing.assert_array_equal(H_0.full(), np.zeros((2, 2), dtype=complex))

    def test_energies_on_diagonal(self):
        """Energies appear on the diagonal of H_0."""
        obs = make_two_state()
        energies = [0.0, 1.5]
        H_0 = obs.set_H_0(energies)
        expected = np.diag(energies).astype(complex)
        np.testing.assert_array_almost_equal(H_0.full(), expected)

    def test_set_H_0_updates_attribute(self):
        """set_H_0 stores the result on self.H_0."""
        obs = make_two_state()
        H_0 = obs.set_H_0([0.0, 2.0])
        self.assertEqual(obs.H_0, H_0)


class TestHIList(unittest.TestCase):
    """Tests for OBBase.H_I_list() and H_I_sum()."""

    def test_H_I_list_returns_H_Omega_list(self):
        """H_I_list() should return H_Omega_list."""
        obs = make_two_state()
        self.assertIs(obs.H_I_list(), obs.H_Omega_list)

    def test_H_I_sum_is_sum_of_H_Omega_list(self):
        """H_I_sum() should equal the sum of all items in H_Omega_list."""
        obs = make_two_state()
        # Replace H_Omega_list with two known matrices.
        A = qu.Qobj(np.array([[0, 1], [1, 0]], dtype=complex))
        B = qu.Qobj(np.array([[1, 0], [0, -1]], dtype=complex))
        obs.H_Omega_list = [A, B]
        total = obs.H_I_sum()
        expected = (A + B).full()
        np.testing.assert_array_almost_equal(total.full(), expected)
