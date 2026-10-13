"""
Unit tests for the module.

Thomas Ogden <t@ogden.eu>

"""

import json
import os
import unittest

import numpy as np

from clerq import PropagationSolver, spectral, t_funcs, utility

# Absolute path of tests/json directory, so that tests can be called from
# different directories.
JSON_DIR = os.path.abspath(os.path.join(__file__, "../", "json"))


class TestInit(unittest.TestCase):
    def test_init_default(self):
        """Default PropagationSolver has a single-state atom, standard time/space grids."""
        mbs = PropagationSolver()
        self.assertEqual(mbs.system.num_states, 1)
        self.assertEqual(mbs.t_min, 0.0)
        self.assertEqual(mbs.t_max, 1.0)
        self.assertEqual(mbs.t_steps, 100)
        self.assertEqual(mbs.z_min, 0.0)
        self.assertEqual(mbs.z_max, 1.0)
        self.assertEqual(mbs.z_steps, 10)

    def test_init_from_json(self):
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        PropagationSolver().from_json(json_path)


class TestPropagationSolver(unittest.TestCase):
    def test_mb_solve(self):
        """Basic test of mb_solve method."""

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mb_solve_00 = PropagationSolver().from_json(json_path)

        mb_solve_00.solve()

    def test_no_atoms(self):
        """Setting the number density ampl to 0.0, i.e. no atoms. The end
        pulse should be the same as the start."""

        json_path = os.path.join(JSON_DIR, "mb_solve_no_atoms.json")

        mbs = PropagationSolver().from_json(json_path)

        mbs.solve(step="euler")

        self.assertEqual(mbs.Omegas_zt.shape, (1, 5, 101))

        # Check that the field at the end of the medium matches the field
        # at the start of the medium.
        self.assertTrue(
            np.allclose(mbs.Omegas_zt[:, 0, :], mbs.Omegas_zt[:, -1, :], rtol=1.0e-6)
        )

    def test_no_atoms_ab(self):
        """Setting the number density to 0.0, i.e. no atoms, with AB step."""

        json_path = os.path.join(JSON_DIR, "mb_solve_no_atoms.json")

        mbs = PropagationSolver().from_json(json_path)

        mbs.solve(step="ab")

        # Check that the field at the end of the medium matches the field
        # at the start of the medium.
        self.assertTrue(
            np.allclose(mbs.Omegas_zt[:, 0, :], mbs.Omegas_zt[:, -1, :], rtol=1.0e-6)
        )

    def test_no_decays(self):
        """Empty decay list."""

        json_path = os.path.join(JSON_DIR, "mb_solve_no_decays.json")

        mb_solve_nd = PropagationSolver().from_json(json_path)

        mb_solve_nd.solve()

    def test_no_rabi_freq_t_func(self):
        """Empty decay list. TODO: No mbsolve, should be in init"""

        json_path = os.path.join(JSON_DIR, "mb_solve_no_rabi_freq_t_func.json")
        mbs = PropagationSolver().from_json(json_path)

        # self.assertEqual(mbs.ob_atom.fields[0].rabi_freq_t_func,
        # t_funcs.square_1)
        self.assertDictEqual(
            mbs.system.fields[0].rabi_freq_t_args,
            {"ampl_0": 1.0, "on_0": 0.0, "off_0": 1.0},
        )

    def test_two_gaussian_2pi(self):
        """Test of a gaussian input 2pi soliton propagating through a two-level
        system.
        """
        json_path = os.path.join(JSON_DIR, "mbs_two_gaussian_2pi.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()
        # Input pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][0] / (np.pi), 2.0, places=1)
        # Output pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][-1] / (np.pi), 2.0, places=1)

    def test_two_gaussian_2pi_n_pi(self):
        """Test of a gaussian input 2pi soliton propagating through a two-level
        system.
        """
        json_path = os.path.join(JSON_DIR, "mbs_two_gaussian_2pi_n_pi.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()
        # Input pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][0] / (np.pi), 2.0, places=1)
        # Output pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][-1] / (np.pi), 2.0, places=1)

    def test_two_sech_2pi(self):
        """Test of a 2pi soliton propagating through a two-level system."""
        json_path = os.path.join(JSON_DIR, "mbs_two_sech_2pi.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()
        # Input pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][0] / (np.pi), 2.0, places=1)
        # Output pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][-1] / (np.pi), 2.0, places=1)

    def test_two_sech_2pi_n_pi(self):
        """Test of a 2pi soliton propagating through a two-level system,
        passing n_pi.
        """
        json_path = os.path.join(JSON_DIR, "mbs_two_sech_2pi_n_pi.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()
        # Input pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][0] / (np.pi), 2.0, places=1)
        # Output pulse is 2pi
        self.assertAlmostEqual(mbs.fields_area()[0][-1] / (np.pi), 2.0, places=1)

    def test_lambda_eit_probe_absorbed_without_coupling(self):
        """Without a coupling field a weak resonant probe is near-completely
        absorbed by the Λ medium (interaction_strength=10, decay rate=1).
        """
        json_path = os.path.join(JSON_DIR, "mbs_lambda_eit_no_coupling.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()

        probe_area = mbs.fields_area()[0]
        transmission = probe_area[-1] / probe_area[0]

        self.assertLess(transmission, 0.05)

    def test_lambda_eit_probe_transmitted_with_coupling(self):
        """With a strong coupling field (Ω_c = 5 >> decay rate = 1) the same
        probe is substantially transmitted through the EIT window.
        """
        json_path = os.path.join(JSON_DIR, "mbs_lambda_eit_with_coupling.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()

        probe_area = mbs.fields_area()[0]
        transmission = probe_area[-1] / probe_area[0]

        self.assertGreater(transmission, 0.4)

    def test_no_vel_classes(self):
        """Empty velocity class dict."""

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)

        vc = {}

        mbs.build_velocity_classes(vc)
        mbs.solve()

    def test_no_vel_classes_inner(self):
        """No inner delta values in dict. TODO: No mbsolve, should be init"""

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)

        vc = {
            "thermal_delta_min": -1.0,
            "thermal_delta_max": 1.0,
            "thermal_delta_steps": 2,
            "thermal_width": 1.0,
        }

        mbs.build_velocity_classes(vc)
        mbs.solve()

    def test_zero_thermal_width(self):
        """TODO: No mbsolve, should be in init"""
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)

        vc = {
            "thermal_delta_min": -1.0,
            "thermal_delta_max": 1.0,
            "thermal_delta_steps": 2,
            "thermal_delta_inner_min": 0.0,
            "thermal_delta_inner_max": 0.0,
            "thermal_delta_inner_steps": 0,
            "thermal_width": 0.0,
        }

        self.assertRaises(ValueError, mbs.build_velocity_classes, vc)

    def test_vel_classes(self):
        """Tests that for a linear two-level system with velocity classes, the
        absorption matches the known Voigt profile.
        """
        json_path = os.path.join(JSON_DIR, "velocity-classes.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()
        freq_list = spectral.freq_list(mbs)
        abs = spectral.absorption(mbs, 0, -1)
        voigt = spectral.voigt_two_linear_known(freq_list, 1.0, 0.05).imag
        # Assert that the max of the abs residuals between the absorption
        # profile and the known broadened Voigt absorption profile for linear
        # two-level systems is below a tolerance. The tolerance scales with
        # OD (doubled after fixing the Maxwell propagation factor-of-2 bug).
        self.assertTrue(np.max(np.abs(abs - voigt)) < 0.1)


class TestVelocityClassUnits(unittest.TestCase):
    """Pins the units of the velocity classes (Doppler broadening).

    ``thermal_delta_*`` and ``thermal_width`` are cyclic frequencies, in the
    same units as ``decay`` rates and Rabi frequencies. ``thermal_delta_list``
    is stored in angular units, and the weights are
    ``exp(-(Delta / (2 pi thermal_width))**2)``, so a velocity class with
    ``thermal_delta = x`` must detune the Hamiltonian by a cyclic ``x``.
    """

    def test_hamiltonian_sees_cyclic_detuning(self):
        """Each class sets H_Delta from a cyclic detuning equal to
        thermal_delta_list / (2 pi), not thermal_delta_list itself."""
        mbs = PropagationSolver().from_json(os.path.join(JSON_DIR, "mb_solve_01.json"))
        mbs.build_velocity_classes(
            {
                "thermal_delta_min": -1.0,
                "thermal_delta_max": 1.0,
                "thermal_delta_steps": 2,
                "thermal_width": 1.0,
            }
        )
        seen = []
        orig = mbs.system.set_H_Delta

        def spy(detunings):
            seen.append(list(detunings))
            return orig(detunings)

        mbs.system.set_H_Delta = spy
        mbs._solve_system = lambda opts=None, save=False: None
        mbs.states_t = lambda: np.zeros(
            (len(mbs.tlist), mbs.system.num_states, mbs.system.num_states),
            dtype=complex,
        )
        mbs._solve_and_average_over_thermal_detunings()

        fixed = mbs.system.get_detunings()[0]
        # One call per velocity class, then one restoring the fixed detunings.
        cyclic = [s[0] - fixed for s in seen[:-1]]
        np.testing.assert_allclose(cyclic, [-1.0, 0.0, 1.0], atol=1e-12)

    def test_absorption_matches_doppler_broadened_lorentzian(self):
        """The weak-probe absorption of a two-level system is the Lorentzian
        (HWHM rate/2) convolved with the Gaussian exp(-(f/w)**2) of the
        cyclic thermal width w. The reference is a direct quadrature, not
        spectral.voigt_two_linear_known."""
        with open(os.path.join(JSON_DIR, "velocity-classes.json")) as fh:
            problem = json.load(fh)
        rate, width, g = 1.0, 0.3, 1.0
        problem["z_steps"] = 10
        problem["interaction_strengths"] = [g]
        problem["velocity_classes"] = {
            "thermal_delta_min": -1.5,
            "thermal_delta_max": 1.5,
            "thermal_delta_steps": 60,
            "thermal_width": width,
        }
        mbs = PropagationSolver().from_json_str(json.dumps(problem))
        mbs.solve()
        freq_list = spectral.freq_list(mbs)
        absorption = spectral.absorption(mbs, 0, -1)

        x = np.linspace(-6.0, 6.0, 4001)
        gauss = np.exp(-((x / width) ** 2)) / (width * np.sqrt(np.pi))
        known = np.array(
            [
                np.trapezoid(
                    gauss * g / (rate / 2) / (1 + ((f - x) / (rate / 2)) ** 2), x
                )
                for f in freq_list
            ]
        )

        self.assertLess(np.max(np.abs(absorption - known)), 0.1)
        # Peak height distinguishes the cyclic width from 2 pi times it.
        self.assertAlmostEqual(absorption.max() / known.max(), 1.0, delta=0.05)
        # The packaged Voigt profile agrees with the quadrature.
        voigt = spectral.voigt_two_linear_known(freq_list, rate, width).imag
        self.assertLess(np.max(np.abs(voigt * g - known)), 1e-3)


class TestSaveLoad(unittest.TestCase):
    """Tests for the PropagationSolver save and load methods."""

    def test_save_load_01(self):
        """Solve a basic PropagationSolver problem. Save the results to file. Set the
        results in the PropagationSolver object to null. Load the results from
        file and check that they equal the original values.
        """
        import tempfile

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mb_solve_01 = PropagationSolver().from_json(json_path)

        Omegas_zt, states_zt = mb_solve_01.solve()

        with tempfile.TemporaryDirectory() as tmpdir:
            mb_solve_01.savefile = os.path.join(tmpdir, "result")
            mb_solve_01.save_results()

            mb_solve_01.Omegas_zt = None
            mb_solve_01.states_zt = None

            mb_solve_01.load_results()

        Omegas_zt_loaded = mb_solve_01.Omegas_zt
        states_zt_loaded = mb_solve_01.states_zt

        self.assertTrue((Omegas_zt == Omegas_zt_loaded).all())
        self.assertTrue((states_zt == states_zt_loaded).all())

    def test_save_load_no_recalc(self):
        import tempfile

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mb_solve_01 = PropagationSolver().from_json(json_path)

        Omegas_zt, states_zt = mb_solve_01.solve()

        with tempfile.TemporaryDirectory() as tmpdir:
            mb_solve_01.savefile = os.path.join(tmpdir, "result")
            mb_solve_01.save_results()

            mb_solve_01.Omegas_zt = None
            mb_solve_01.states_zt = None

            Omegas_zt, states_zt = mb_solve_01.solve(recalc=False)

        Omegas_zt_loaded = mb_solve_01.Omegas_zt
        states_zt_loaded = mb_solve_01.states_zt

        self.assertTrue((Omegas_zt == Omegas_zt_loaded).all())
        self.assertTrue((states_zt == states_zt_loaded).all())

    def test_save_creates_missing_directory(self):
        """GH#195: save_results must not error when the savefile directory
        does not yet exist."""
        import tempfile

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()

        with tempfile.TemporaryDirectory() as tmpdir:
            savefile = os.path.join(tmpdir, "subdir", "nested", "result")
            mbs.savefile = savefile
            mbs.save_results()  # must not raise
            self.assertTrue(os.path.isfile(savefile + ".qu"))

    def test_load_hash_mismatch_raises(self):
        """GH#199: load_results raises ValueError when problem definition
        has changed since the savefile was written."""
        import tempfile

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs_a = PropagationSolver().from_json(json_path)
        mbs_a.solve()

        # Build a different problem (more z_steps) and reuse the savefile from mbs_a.
        mbs_b = PropagationSolver().from_json(json_path)
        mbs_b.z_steps = mbs_a.z_steps + 1
        mbs_b.build_zlist(mbs_b.z_min, mbs_b.z_max, mbs_b.z_steps, mbs_b.z_steps_inner)

        with tempfile.TemporaryDirectory() as tmpdir:
            savefile = os.path.join(tmpdir, "result")
            mbs_a.savefile = savefile
            mbs_a.save_results()

            mbs_b.savefile = savefile
            with self.assertRaises(ValueError):
                mbs_b.load_results()

    def test_load_stale_file_resolves(self):
        """GH#199: mbsolve() re-solves automatically when the savefile hash
        does not match the current problem definition."""
        import tempfile

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs_a = PropagationSolver().from_json(json_path)
        mbs_a.solve()

        # Build a different problem and point it at mbs_a's stale savefile.
        mbs_b = PropagationSolver().from_json(json_path)
        mbs_b.z_steps = mbs_a.z_steps + 1
        mbs_b.build_zlist(mbs_b.z_min, mbs_b.z_max, mbs_b.z_steps, mbs_b.z_steps_inner)

        with tempfile.TemporaryDirectory() as tmpdir:
            savefile = os.path.join(tmpdir, "result")
            mbs_a.savefile = savefile
            mbs_a.save_results()

            mbs_b.savefile = savefile
            Omegas_zt, _ = mbs_b.solve(recalc=False)

        # Shape should reflect mbs_b's z_steps, not the stale mbs_a result.
        self.assertEqual(Omegas_zt.shape[1], mbs_b.z_steps + 1)

    def test_load_old_format_warns(self):
        """GH#199: loading a legacy 2-tuple savefile emits a warning but
        still loads the results."""
        import tempfile
        import warnings

        import qutip as qu

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.solve()

        with tempfile.TemporaryDirectory() as tmpdir:
            savefile = os.path.join(tmpdir, "result")
            mbs.savefile = savefile
            # Write old-format 2-tuple directly.
            qu.qsave((mbs.Omegas_zt, mbs.states_zt), savefile)

            mbs.Omegas_zt = None
            mbs.states_zt = None

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                mbs.load_results()

        self.assertTrue(any("no integrity metadata" in str(w.message) for w in caught))
        self.assertIsNotNone(mbs.Omegas_zt)


class TestBuildZlist(unittest.TestCase):
    def test_default_zlist(self):
        mbs = PropagationSolver()
        zlist = np.linspace(0.0, 1.0, 11)
        self.assertTrue(np.allclose(mbs.zlist, zlist, rtol=1.0e-6))

    def test_z_step(self):
        """z_step returns (z_max - z_min) / z_steps."""
        mbs = PropagationSolver()
        self.assertAlmostEqual(mbs.z_step(), 0.1)

    def test_z_step_inner(self):
        """z_step_inner divides z_step by z_steps_inner."""
        mbs = PropagationSolver()
        self.assertAlmostEqual(mbs.z_step_inner(), mbs.z_step() / mbs.z_steps_inner)


class TestCheck(unittest.TestCase):
    def test_check_passes_when_strengths_match_fields(self):
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        self.assertTrue(mbs.check())

    def test_check_raises_when_strengths_mismatch_fields(self):
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        mbs.interaction_strengths = [1.0, 2.0]  # two strengths, one field
        self.assertRaises(ValueError, mbs.check)


class TestBuildIntpHOmegaList(unittest.TestCase):
    """Unit tests for PropagationSolver._build_intp_H_Omega_list."""

    def test_one_field_structure(self):
        """Returns a list of [Qobj, Coefficient] pairs, one per field."""
        import qutip as qu

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        Omegas_z = mbs.Omegas_zt[:, 0, :]
        result = mbs._build_intp_H_Omega_list(Omegas_z)

        self.assertEqual(len(result), 1)
        H_op, coeff = result[0]
        from qutip.core.cy.coefficient import Coefficient

        self.assertIsInstance(H_op, qu.Qobj)
        self.assertIsInstance(coeff, Coefficient)

    def test_two_fields_length(self):
        """Returns one pair per field."""
        json_path = os.path.join(JSON_DIR, "mb_solve_lamda.json")
        mbs = PropagationSolver().from_json(json_path)
        Omegas_z = mbs.Omegas_zt[:, 0, :]
        result = mbs._build_intp_H_Omega_list(Omegas_z)

        self.assertEqual(len(result), 2)

    def test_coefficient_evaluates_correctly(self):
        """Coefficient at tlist[k] equals Omegas_z[0,k].real / (2π)."""
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        # Use a non-trivial field profile
        mbs.solve()
        Omegas_z = mbs.Omegas_zt[:, -1, :]
        result = mbs._build_intp_H_Omega_list(Omegas_z)
        _H_op, coeff = result[0]
        # Check a few time points
        for k in [0, len(mbs.tlist) // 2, -1]:
            t = mbs.tlist[k]
            expected = Omegas_z[0, k].real / (2 * np.pi)
            self.assertAlmostEqual(coeff(t), expected, places=10)


class TestPopulations(unittest.TestCase):
    def test_twolevel_shape(self):

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)

        pop_lower = mbs.populations([0])
        pop_upper = mbs.populations([1])

        np.testing.assert_allclose(
            pop_lower, np.zeros((mbs.z_steps + 1, mbs.t_steps + 1))
        )
        np.testing.assert_allclose(
            pop_upper, np.zeros((mbs.z_steps + 1, mbs.t_steps + 1))
        )


class TestPopulationsField(unittest.TestCase):
    def test_twolevel_shape(self):

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)

        pop_upper = mbs.populations_field(field_idx=0, upper=True)
        pop_lower = mbs.populations_field(field_idx=0, upper=False)

        np.testing.assert_allclose(
            pop_lower, np.zeros((mbs.z_steps + 1, mbs.t_steps + 1))
        )
        np.testing.assert_allclose(
            pop_upper, np.zeros((mbs.z_steps + 1, mbs.t_steps + 1))
        )


class TestCoherences(unittest.TestCase):
    def test_twolevel_shape(self):

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        coh = mbs.coherences([[0, 1]])

        np.testing.assert_allclose(coh, np.zeros((mbs.z_steps + 1, mbs.t_steps + 1)))


class TestCoherencesField(unittest.TestCase):
    def test_twolevel_shape(self):

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        coh = mbs.coherences_field(field_idx=0)

        np.testing.assert_allclose(coh, np.zeros((mbs.z_steps + 1, mbs.t_steps + 1)))


class TestNormaliseVelocityClasses(unittest.TestCase):
    """Unit tests for PropagationSolver._normalise_velocity_classes."""

    def setUp(self):
        self.mbs = PropagationSolver()

    def test_none_gives_all_defaults(self):
        """None input returns a dict with all default keys."""
        vc = self.mbs._normalise_velocity_classes(None)
        self.assertEqual(vc["thermal_width"], 1.0)
        self.assertEqual(vc["thermal_delta_min"], 0.0)
        self.assertEqual(vc["thermal_delta_max"], 0.0)
        self.assertEqual(vc["thermal_delta_steps"], 0)
        self.assertEqual(vc["thermal_delta_inner_min"], 0.0)
        self.assertEqual(vc["thermal_delta_inner_max"], 0.0)
        self.assertEqual(vc["thermal_delta_inner_steps"], 0)

    def test_empty_dict_gives_all_defaults(self):
        """Empty dict returns a dict with all default keys."""
        vc = self.mbs._normalise_velocity_classes({})
        self.assertEqual(vc["thermal_width"], 1.0)
        self.assertEqual(vc["thermal_delta_steps"], 0)

    def test_user_values_override_defaults(self):
        """Supplied values override the defaults; unspecified keys keep defaults."""
        vc = self.mbs._normalise_velocity_classes(
            {"thermal_width": 2.5, "thermal_delta_min": -1.0, "thermal_delta_max": 1.0}
        )
        self.assertEqual(vc["thermal_width"], 2.5)
        self.assertEqual(vc["thermal_delta_min"], -1.0)
        self.assertEqual(vc["thermal_delta_max"], 1.0)
        # Unspecified keys still present with defaults
        self.assertEqual(vc["thermal_delta_steps"], 0)
        self.assertEqual(vc["thermal_delta_inner_steps"], 0)

    def test_zero_thermal_width_raises(self):
        with self.assertRaises(ValueError):
            self.mbs._normalise_velocity_classes({"thermal_width": 0.0})

    def test_negative_thermal_width_raises(self):
        with self.assertRaises(ValueError):
            self.mbs._normalise_velocity_classes({"thermal_width": -1.0})

    def test_does_not_mutate_input(self):
        """The original dict must not be modified."""
        original = {"thermal_width": 0.5}
        original_copy = dict(original)
        self.mbs._normalise_velocity_classes(original)
        self.assertEqual(original, original_copy)


class TestBuildThermalDeltaList(unittest.TestCase):
    """Unit tests for PropagationSolver._build_thermal_delta_list."""

    def setUp(self):
        self.mbs = PropagationSolver()

    def _vc(self, **kwargs):
        """Helper: build a normalised vc dict from kwargs."""
        return self.mbs._normalise_velocity_classes(kwargs if kwargs else None)

    def test_single_point_default(self):
        """Default vc (all zeros) gives a single detuning at 0."""
        vc = self._vc()
        deltas = self.mbs._build_thermal_delta_list(vc)
        np.testing.assert_array_almost_equal(deltas, [0.0])

    def test_outer_range_only(self):
        """Outer range with 2 steps gives 3 evenly-spaced points."""
        vc = self._vc(
            thermal_delta_min=-1.0,
            thermal_delta_max=1.0,
            thermal_delta_steps=2,
        )
        deltas = self.mbs._build_thermal_delta_list(vc)
        expected = 2 * np.pi * np.array([-1.0, 0.0, 1.0])
        np.testing.assert_array_almost_equal(deltas, expected)

    def test_two_pi_scaling(self):
        """Detuning values are in angular-frequency units (multiplied by 2π)."""
        vc = self._vc(
            thermal_delta_min=0.0,
            thermal_delta_max=1.0,
            thermal_delta_steps=1,
        )
        deltas = self.mbs._build_thermal_delta_list(vc)
        np.testing.assert_array_almost_equal(deltas, [0.0, 2 * np.pi])

    def test_inner_and_outer_merged_and_deduplicated(self):
        """Overlapping outer and inner ranges are merged with np.unique."""
        vc = self._vc(
            thermal_delta_min=-1.0,
            thermal_delta_max=1.0,
            thermal_delta_steps=2,
            thermal_delta_inner_min=-0.5,
            thermal_delta_inner_max=0.5,
            thermal_delta_inner_steps=2,
        )
        deltas = self.mbs._build_thermal_delta_list(vc)
        # Outer: -1, 0, 1  Inner: -0.5, 0, 0.5  — 0 appears in both, deduped
        expected = 2 * np.pi * np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
        np.testing.assert_array_almost_equal(deltas, expected)

    def test_output_is_sorted(self):
        """Output array must be monotonically increasing."""
        vc = self._vc(
            thermal_delta_min=-2.0,
            thermal_delta_max=2.0,
            thermal_delta_steps=4,
            thermal_delta_inner_min=-0.5,
            thermal_delta_inner_max=0.5,
            thermal_delta_inner_steps=2,
        )
        deltas = self.mbs._build_thermal_delta_list(vc)
        self.assertTrue(np.all(np.diff(deltas) > 0))


class TestZStepFields(unittest.TestCase):
    """Unit tests for _z_step_fields_euler and _z_step_fields_ab."""

    def setUp(self):
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        self.mbs = PropagationSolver().from_json(json_path)

    def test_euler_output_shape(self):
        n_fields = len(self.mbs.system.fields)
        n_t = len(self.mbs.tlist)
        Omegas_z = np.zeros((n_fields, n_t), dtype=complex)
        sum_coh = np.zeros((n_fields, n_t), dtype=complex)
        result = self.mbs._z_step_fields_euler(
            h=0.1, N=1.0, Omegas_z_this=Omegas_z, sum_coh_this=sum_coh
        )
        self.assertEqual(result.shape, (n_fields, n_t))

    def test_euler_zero_coherence_is_identity(self):
        """With zero coherence the field does not change."""
        n_fields = len(self.mbs.system.fields)
        n_t = len(self.mbs.tlist)
        Omegas_z = np.ones((n_fields, n_t), dtype=complex)
        sum_coh = np.zeros((n_fields, n_t), dtype=complex)
        result = self.mbs._z_step_fields_euler(
            h=0.1, N=1.0, Omegas_z_this=Omegas_z, sum_coh_this=sum_coh
        )
        np.testing.assert_array_equal(result, Omegas_z)

    def test_euler_known_value(self):
        """dΩ/dz = 2i·N·g·ρ; check against hand-computed value."""
        n_fields = len(self.mbs.system.fields)
        n_t = len(self.mbs.tlist)
        Omegas_z = np.zeros((n_fields, n_t), dtype=complex)
        sum_coh = np.ones((n_fields, n_t), dtype=complex)
        h, N = 0.1, 1.0
        result = self.mbs._z_step_fields_euler(
            h=h, N=N, Omegas_z_this=Omegas_z, sum_coh_this=sum_coh
        )
        expected = h * 2.0j * N * self.mbs.g[0] * np.ones(n_t, dtype=complex)
        np.testing.assert_allclose(result[0], expected)

    def test_ab_output_shape(self):
        n_fields = len(self.mbs.system.fields)
        n_t = len(self.mbs.tlist)
        Omegas_z = np.zeros((n_fields, n_t), dtype=complex)
        sum_coh = np.zeros((n_fields, n_t), dtype=complex)
        result = self.mbs._z_step_fields_ab(
            h=0.1,
            N=1.0,
            sum_coh_prev=sum_coh,
            sum_coh_this=sum_coh,
            Omegas_z_this=Omegas_z,
        )
        self.assertEqual(result.shape, (n_fields, n_t))

    def test_ab_has_no_z_prev(self):
        import inspect

        sig = inspect.signature(self.mbs._z_step_fields_ab)
        self.assertNotIn("z_prev", sig.parameters)
