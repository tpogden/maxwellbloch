"""
Unit tests for the module.

Thomas Ogden <t@ogden.eu>

"""

import json
import os
import tempfile
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

        self.assertEqual(mbs._Omegas_zt.shape, (1, 5, 101))

        # Check that the field at the end of the medium matches the field
        # at the start of the medium.
        self.assertTrue(
            np.allclose(mbs._Omegas_zt[:, 0, :], mbs._Omegas_zt[:, -1, :], rtol=1.0e-6)
        )

    def test_no_atoms_ab(self):
        """Setting the number density to 0.0, i.e. no atoms, with AB step."""

        json_path = os.path.join(JSON_DIR, "mb_solve_no_atoms.json")

        mbs = PropagationSolver().from_json(json_path)

        mbs.solve(step="ab")

        # Check that the field at the end of the medium matches the field
        # at the start of the medium.
        self.assertTrue(
            np.allclose(mbs._Omegas_zt[:, 0, :], mbs._Omegas_zt[:, -1, :], rtol=1.0e-6)
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

            mb_solve_01._Omegas_zt = None
            mb_solve_01.states_zt = None

            mb_solve_01.load_results()

        Omegas_zt_loaded = mb_solve_01._Omegas_zt
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

            mb_solve_01._Omegas_zt = None
            mb_solve_01.states_zt = None

            Omegas_zt, states_zt = mb_solve_01.solve(recalc=False)

        Omegas_zt_loaded = mb_solve_01._Omegas_zt
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
            qu.qsave((mbs._Omegas_zt, mbs.states_zt), savefile)

            mbs._Omegas_zt = None
            mbs.states_zt = None

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                mbs.load_results()

        self.assertTrue(any("no integrity metadata" in str(w.message) for w in caught))
        self.assertIsNotNone(mbs._Omegas_zt)


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
        """Returns [Qobj, Coefficient] pairs: a real field needs no Im term."""
        import qutip as qu

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        Omegas_z = mbs._Omegas_zt[:, 0, :]
        result = mbs._build_intp_H_Omega_list(Omegas_z)

        self.assertEqual(len(result), 1)
        from qutip.core.cy.coefficient import Coefficient

        for H_op, coeff in result:
            self.assertIsInstance(H_op, qu.Qobj)
            self.assertIsInstance(coeff, Coefficient)

    def test_complex_field_adds_an_imaginary_term(self):
        """A field with an imaginary part gets a second term; a real one does not."""
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver.from_json(json_path)
        Omegas_z = mbs._Omegas_zt[:, 0, :] * (1 + 0.5j)
        result = mbs._build_intp_H_Omega_list(Omegas_z)
        self.assertEqual(len(result), 2)
        _, coeff_im = result[1]
        k = len(mbs.tlist) // 2
        self.assertAlmostEqual(
            coeff_im(mbs.tlist[k]), Omegas_z[0, k].imag / (2 * np.pi), places=10
        )

    def test_real_drive_has_one_term_per_field(self):
        """The deprecated real_drive=True keeps one Re(Omega) term per field."""
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver.from_json(json_path)
        mbs.real_drive = True
        result = mbs._build_intp_H_Omega_list(mbs._Omegas_zt[:, 0, :])
        self.assertEqual(len(result), 1)

    def test_two_fields_length(self):
        """Returns one pair per field for real fields."""
        json_path = os.path.join(JSON_DIR, "mb_solve_lamda.json")
        mbs = PropagationSolver().from_json(json_path)
        Omegas_z = mbs._Omegas_zt[:, 0, :]
        result = mbs._build_intp_H_Omega_list(Omegas_z)

        self.assertEqual(len(result), 2)

    def test_coefficient_evaluates_correctly(self):
        """Coefficient at tlist[k] equals Omegas_z[0,k].real / (2π)."""
        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mbs = PropagationSolver().from_json(json_path)
        # Use a non-trivial field profile
        mbs.solve()
        Omegas_z = mbs._Omegas_zt[:, -1, :]
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


class TestThermalWeights(unittest.TestCase):
    """thermal_weights is the Boltzmann density; the Doppler average uses it
    times the detuning interval each class represents (thermal_widths)."""

    BANDED = {
        "thermal_width": 1.0,
        "thermal_delta_min": -4.0,
        "thermal_delta_max": 4.0,
        "thermal_delta_steps": 8,
        "thermal_delta_inner_min": -1.0,
        "thermal_delta_inner_max": 1.0,
        "thermal_delta_inner_steps": 6,
    }

    def setUp(self):
        self.mbs = PropagationSolver()

    def _quadrature(self, **vc):
        self.mbs.build_velocity_classes({"thermal_width": 1.0, **vc})
        q = self.mbs.thermal_quadrature_weights
        return self.mbs.thermal_delta_list, q / q.sum()

    def test_thermal_weights_stay_the_boltzmann_density(self):
        deltas, w = self.mbs.build_velocity_classes(self.BANDED)
        np.testing.assert_allclose(
            w, utility.maxwell_boltzmann(deltas, 2 * np.pi), rtol=1e-12
        )

    def test_uniform_grid_quadrature_is_proportional_to_density(self):
        _, w = self.mbs.build_velocity_classes(
            {
                "thermal_width": 1.0,
                "thermal_delta_min": -3.0,
                "thermal_delta_max": 3.0,
                "thermal_delta_steps": 12,
            }
        )
        q = self.mbs.thermal_quadrature_weights
        np.testing.assert_allclose(q / q.sum(), w / w.sum(), rtol=1e-12)

    def test_widths(self):
        """Interior widths are half the distance between neighbours; the end
        classes take their whole adjacent gap."""
        np.testing.assert_allclose(
            PropagationSolver._thermal_spacings(np.array([0.0, 1.0, 3.0, 4.0])),
            [1.0, 1.5, 1.5, 1.0],
        )

    def test_single_class_quadrature_weight_is_the_density(self):
        deltas, w = self.mbs.build_velocity_classes(None)
        self.assertEqual(len(w), 1)
        np.testing.assert_allclose(self.mbs.thermal_widths, [1.0])
        np.testing.assert_allclose(
            self.mbs.thermal_quadrature_weights,
            utility.maxwell_boltzmann(deltas, 2 * np.pi),
        )

    def test_nonuniform_grid_matches_fine_uniform_grid(self):
        """An inner band must not over-weight the dense region."""
        fine = {
            "thermal_delta_min": -6.0,
            "thermal_delta_max": 6.0,
            "thermal_delta_steps": 2400,
        }
        banded = {
            "thermal_delta_min": -6.0,
            "thermal_delta_max": 6.0,
            "thermal_delta_steps": 120,
            "thermal_delta_inner_min": -1.0,
            "thermal_delta_inner_max": 1.0,
            "thermal_delta_inner_steps": 40,
        }
        d_f, w_f = self._quadrature(**fine)
        d_b, w_b = self._quadrature(**banded)
        # Compare a smooth observable: the weighted second moment of Delta.
        m2_f = np.sum(w_f * d_f**2)
        m2_b = np.sum(w_b * d_b**2)
        self.assertLess(abs(m2_b / m2_f - 1.0), 2e-3)
        # The unweighted Boltzmann density over-weights the dense inner band.
        raw = utility.maxwell_boltzmann(d_b, 2 * np.pi)
        m2_raw = np.sum(raw / raw.sum() * d_b**2)
        self.assertGreater(abs(m2_raw / m2_f - 1.0), 10 * abs(m2_b / m2_f - 1.0))

    def test_symmetric_grid_gives_symmetric_weights(self):
        _, w = self._quadrature(**self.BANDED)
        np.testing.assert_allclose(w, w[::-1], rtol=1e-9)

    def test_symmetric_banded_grid_keeps_drives_equivalent(self):
        """A mirror-symmetric banded grid has mirrored quadrature weights, so
        _drive_can_differ must not flip because of them."""
        mbs = PropagationSolver().from_json(os.path.join(JSON_DIR, "mb_solve_01.json"))
        mbs.build_velocity_classes(self.BANDED)
        self.assertFalse(mbs._drive_can_differ())
        mbs.build_velocity_classes({**self.BANDED, "thermal_delta_max": 4.5})
        self.assertTrue(mbs._drive_can_differ())

    def _banded_savefile(self, tmpdir, banded, marker):
        import qutip as qu

        mbs = PropagationSolver().from_json(os.path.join(JSON_DIR, "mb_solve_01.json"))
        mbs.build_velocity_classes(self.BANDED if banded else None)
        mbs.solve()
        meta = mbs._savefile_metadata()
        if not marker:
            del meta["thermal_quadrature"]
        mbs.savefile = os.path.join(tmpdir, "result")
        qu.qsave((mbs._Omegas_zt, mbs.states_zt, meta), mbs.savefile)
        return mbs

    def test_savefile_without_marker_rejected_on_nonuniform_grid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mbs = self._banded_savefile(tmpdir, banded=True, marker=False)
            with self.assertRaisesRegex(ValueError, "ignored the class spacing"):
                mbs.load_results()

    def test_savefile_with_marker_loads_on_nonuniform_grid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mbs = self._banded_savefile(tmpdir, banded=True, marker=True)
            mbs.load_results()

    def test_savefile_without_marker_loads_on_uniform_grid(self):
        """Weights are unchanged on a uniform grid, so old results stand."""
        with tempfile.TemporaryDirectory() as tmpdir:
            mbs = self._banded_savefile(tmpdir, banded=False, marker=False)
            mbs.load_results()

    def test_banded_absorption_matches_doppler_broadened_lorentzian(self):
        """Solve-level check: a coarse outer grid plus a fine inner band
        reproduces the Voigt absorption of a two-level system."""
        with open(os.path.join(JSON_DIR, "velocity-classes.json")) as fh:
            problem = json.load(fh)
        rate, width, g = 1.0, 0.3, 1.0
        problem["z_steps"] = 10
        problem["interaction_strengths"] = [g]
        problem["velocity_classes"] = {
            "thermal_delta_min": -1.5,
            "thermal_delta_max": 1.5,
            "thermal_delta_steps": 15,
            "thermal_delta_inner_min": -0.5,
            "thermal_delta_inner_max": 0.5,
            "thermal_delta_inner_steps": 40,
            "thermal_width": width,
        }
        mbs = PropagationSolver().from_json_str(json.dumps(problem))
        mbs.solve()
        freq_list = spectral.freq_list(mbs)
        absorption = spectral.absorption(mbs, 0, -1)
        voigt = spectral.voigt_two_linear_known(freq_list, rate, width).imag * g
        self.assertLess(np.max(np.abs(absorption - voigt)), 0.1)


class TestThermalDeltaSpacing(unittest.TestCase):
    """thermal_delta_spacing builds a uniform grid free of Doppler echoes."""

    def _solver(self, field_extra=(), doppler=None, rate=None, **vc):
        """Two-level problem; ``doppler`` sets factor_doppler_shift of the
        field, and each dict in ``field_extra`` adds another field."""
        with open(os.path.join(JSON_DIR, "velocity-classes.json")) as fh:
            problem = json.load(fh)
        problem["z_steps"] = 4
        fields = problem["system"]["fields"]
        if rate is not None:
            problem["system"]["decays"][0]["rate"] = rate
        if doppler is not None:
            fields[0]["factor_doppler_shift"] = doppler
        for extra in field_extra:
            fields.append({**fields[0], **extra})
        problem["velocity_classes"] = {
            "thermal_width": 0.3,
            "thermal_delta_min": -1.5,
            "thermal_delta_max": 1.5,
            **vc,
        }
        return PropagationSolver().from_json_str(json.dumps(problem))

    def test_float_spacing_gives_a_uniform_grid_no_wider_than_requested(self):
        mbs = self._solver(thermal_delta_spacing=0.07)
        cyclic = mbs.thermal_delta_list / (2 * np.pi)
        self.assertAlmostEqual(cyclic[0], -1.5)
        self.assertAlmostEqual(cyclic[-1], 1.5)
        gaps = np.diff(cyclic)
        np.testing.assert_allclose(gaps, gaps[0])
        self.assertLessEqual(gaps[0], 0.07)
        self.assertGreater(gaps[0], 0.07 * 0.97)

    def test_exact_multiple_does_not_add_a_class(self):
        mbs = self._solver(thermal_delta_spacing=0.25)
        self.assertEqual(len(mbs.thermal_delta_list), 13)

    def test_auto_spacing_is_set_by_the_time_window(self):
        """With a slow decay the echo is not damped, so the window sets it."""
        mbs = self._solver(rate=0.05, thermal_delta_spacing="auto")
        window = mbs.t_max - mbs.t_min
        expected = 0.75 / window
        self.assertAlmostEqual(mbs.thermal_delta_auto_spacing(), expected)
        gaps = np.diff(mbs.thermal_delta_list) / (2 * np.pi)
        self.assertLessEqual(gaps.max(), expected * (1 + 1e-9))
        self.assertGreater(gaps.max(), expected * 0.97)

    def test_auto_spacing_shrinks_with_the_doppler_slope(self):
        """A coherence oscillating at b * Delta echoes at 1 / (b h), so the
        spacing divides by b, the difference of the levels' Doppler slopes.
        Fields on one transition add, so a counter-propagating twin cancels."""
        base = self._solver(rate=0.05, thermal_delta_spacing="auto")
        window = base.t_max - base.t_min
        cases = [
            (dict(doppler=2.5), 2.5),
            (dict(doppler=-1.0), 1.0),
            (dict(field_extra=[{"factor_doppler_shift": 1.0}]), 2.0),
            (dict(doppler=0.5, field_extra=[{"factor_doppler_shift": 2.0}]), 2.5),
            # slopes cancel: no Doppler dependence, so the window rule with b = 1
            (dict(field_extra=[{"factor_doppler_shift": -1.0}]), 1.0),
        ]
        for kwargs, slope in cases:
            with self.subTest(**{k: str(v) for k, v in kwargs.items()}):
                mbs = self._solver(rate=0.05, thermal_delta_spacing="auto", **kwargs)
                self.assertAlmostEqual(
                    mbs.thermal_delta_auto_spacing(), 0.75 / (window * slope)
                )

    def test_auto_spacing_relaxes_with_the_decay(self):
        """A fast decay damps the echo below the tolerance, so the spacing is
        set by the decay, g / (b ln(1 / tol)), once that exceeds the window's."""
        mbs = self._solver(rate=1.0, thermal_delta_spacing="auto")
        window = mbs.t_max - mbs.t_min
        by_decay = np.pi / np.log(1e6)  # g = pi * rate for a two-level decay
        self.assertGreater(by_decay, 0.75 / window)
        self.assertAlmostEqual(mbs.thermal_delta_auto_spacing(), by_decay)
        self.assertAlmostEqual(
            mbs.thermal_delta_auto_spacing(tolerance=1e-3), np.pi / np.log(1e3)
        )

    def test_echo_tolerance_key(self):
        mbs = self._solver(
            rate=1.0, thermal_delta_spacing="auto", thermal_delta_echo_tolerance=1e-3
        )
        gaps = np.diff(mbs.thermal_delta_list) / (2 * np.pi)
        self.assertLessEqual(gaps.max(), np.pi / np.log(1e3) * (1 + 1e-9))
        self.assertGreater(gaps.max(), np.pi / np.log(1e3) * 0.9)
        for bad in (0.0, 1.0, 2.0, -1e-3, True, "x"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self._solver(
                    thermal_delta_spacing="auto", thermal_delta_echo_tolerance=bad
                )
        with self.assertRaises(ValueError):  # only with "auto"
            self._solver(thermal_delta_spacing=0.1, thermal_delta_echo_tolerance=1e-3)

    def test_ladder_level_slopes_decays_and_spacing(self):
        """Fields add their detuning to their own upper level, so a ladder's
        levels move at -s_A and -s_B, and the (1, 2) coherence at their
        difference."""

        def ladder(factor_b):
            system = {
                "num_states": 3,
                "decays": [
                    {"channels": [[0, 1]], "rate": 1.0},
                    {"channels": [[1, 2]], "rate": 0.001},
                ],
                "fields": [
                    {"coupled_levels": [[0, 1]], "rabi_freq": 1e-3},
                    {
                        "coupled_levels": [[1, 2]],
                        "rabi_freq": 0.0,
                        "factor_doppler_shift": factor_b,
                    },
                ],
            }
            return PropagationSolver(
                system=system,
                t_min=-1.0,
                t_max=12.0,
                t_steps=100,
                velocity_classes={
                    "thermal_delta_min": -9.0,
                    "thermal_delta_max": 9.0,
                    "thermal_delta_spacing": "auto",
                },
            )

        co = ladder(1.0)
        np.testing.assert_allclose(co._level_doppler_slopes(), [0.0, -1.0, -1.0])
        np.testing.assert_allclose(
            co._level_decay_rates(), [0.0, np.pi * 1.0, np.pi * 0.001]
        )
        # Pairs (0, 1) and (0, 2) have b = 1; (0, 2) is nearly undamped, so the
        # window rule governs; (1, 2) has b = 0 and is skipped.
        self.assertAlmostEqual(co.thermal_delta_auto_spacing(), 0.75 / 13.0)

        counter = ladder(-1.625)
        np.testing.assert_allclose(counter._level_doppler_slopes(), [0.0, -1.0, 1.625])
        # (0, 2) has b = 1.625 and no damping: 0.75 / (13 * 1.625), the smallest.
        self.assertAlmostEqual(
            counter.thermal_delta_auto_spacing(), 0.75 / (13 * 1.625)
        )

    def test_auto_grid_is_converged_with_a_doppler_factor(self):
        """With factor 2 the echo is at 1 / (2 h); the bare 0.75 / window
        spacing puts it inside the window and fails, 'auto' does not."""
        # A slow decay, so that the echo is not damped away.
        fine = self._solver(doppler=2.0, rate=0.05, thermal_delta_steps=400)
        fine.solve()
        ref = np.asarray(fine._Omegas_zt[-1])
        scale = np.abs(ref).max()

        auto = self._solver(doppler=2.0, rate=0.05, thermal_delta_spacing="auto")
        auto.solve()
        err_auto = np.abs(np.asarray(auto._Omegas_zt[-1]) - ref).max() / scale
        self.assertLess(err_auto, 1e-5)

        bare = self._solver(
            doppler=2.0,
            rate=0.05,
            thermal_delta_spacing=0.75 / (auto.t_max - auto.t_min),
        )
        bare.solve()
        err_bare = np.abs(np.asarray(bare._Omegas_zt[-1]) - ref).max() / scale
        self.assertGreater(err_bare, 100 * err_auto)

    def test_numpy_scalars_are_accepted_and_bool_is_not(self):
        mbs = self._solver(thermal_delta_spacing=0.1)
        for good in (np.float32(0.1), np.float64(0.1), np.int64(1)):
            with self.subTest(good=repr(good)):
                mbs.build_velocity_classes(
                    {
                        **mbs.velocity_classes,
                        "thermal_delta_spacing": good,
                    }
                )
                self.assertIsInstance(
                    mbs.velocity_classes["thermal_delta_spacing"], float
                )
                mbs.to_json_str()  # still JSON-serialisable for the hash
        with self.assertRaises(ValueError):
            self._solver(thermal_delta_spacing=True)

    def test_tiny_spacing_warns_about_the_class_count(self):
        with self.assertWarnsRegex(UserWarning, "velocity classes"):
            self._solver(thermal_delta_spacing=1e-3)

    def test_auto_survives_a_json_round_trip(self):
        mbs = self._solver(thermal_delta_spacing="auto")
        again = PropagationSolver().from_json_str(mbs.to_json_str())
        self.assertEqual(again.velocity_classes["thermal_delta_spacing"], "auto")
        np.testing.assert_allclose(again.thermal_delta_list, mbs.thermal_delta_list)

    def test_unused_key_leaves_the_problem_hash_unchanged(self):
        mbs = self._solver(thermal_delta_steps=10)
        self.assertNotIn("thermal_delta_spacing", mbs.velocity_classes)

    def test_spacing_with_steps_raises(self):
        for key in ("thermal_delta_steps", "thermal_delta_inner_steps"):
            with self.subTest(key=key), self.assertRaises(ValueError):
                self._solver(thermal_delta_spacing=0.1, **{key: 4})

    def test_bad_spacing_raises(self):
        for bad in (0.0, -0.1, "fine"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self._solver(thermal_delta_spacing=bad)

    def test_empty_range_raises(self):
        with self.assertRaises(ValueError):
            self._solver(
                thermal_delta_spacing="auto",
                thermal_delta_min=0.0,
                thermal_delta_max=0.0,
            )

    def test_auto_grid_is_converged_and_a_coarser_one_is_not(self):
        """With a slow decay, the echo-free grid matches a fine reference and a
        spacing of 2 / window, with its echo inside the window, fails."""
        fine = self._solver(rate=0.05, thermal_delta_steps=240)
        fine.solve()
        ref = np.asarray(fine._Omegas_zt[-1])
        scale = np.abs(ref).max()

        auto = self._solver(rate=0.05, thermal_delta_spacing="auto")
        auto.solve()
        err_auto = np.abs(np.asarray(auto._Omegas_zt[-1]) - ref).max() / scale
        self.assertLess(err_auto, 1e-4)

        window = auto.t_max - auto.t_min
        coarse = self._solver(rate=0.05, thermal_delta_spacing=2.0 / window)
        coarse.solve()
        err_coarse = np.abs(np.asarray(coarse._Omegas_zt[-1]) - ref).max() / scale
        self.assertGreater(err_coarse, 100 * err_auto)

    def test_decay_limited_auto_grid_is_converged_with_fewer_classes(self):
        """A fast decay damps the echo, so 'auto' uses a coarser grid than the
        window rule and is still converged; ignoring the decay would be wrong
        only for the slow case above."""
        fine = self._solver(rate=1.0, thermal_delta_steps=240)
        fine.solve()
        ref = np.asarray(fine._Omegas_zt[-1])
        scale = np.abs(ref).max()

        auto = self._solver(rate=1.0, thermal_delta_spacing="auto")
        window = auto.t_max - auto.t_min
        self.assertGreater(
            np.diff(auto.thermal_delta_list).min() / (2 * np.pi), 2 * 0.75 / window
        )
        auto.solve()
        err = np.abs(np.asarray(auto._Omegas_zt[-1]) - ref).max() / scale
        self.assertLess(err, 1e-4)


class TestInnerSteps(unittest.TestCase):
    """z_steps_inner subdivides each recorded z step; it must not shift slices.

    Regression: with step="ab" the Euler bootstrap used to store slice 1 after a
    single inner step, so every slice was one inner step (h_inner) early whenever
    z_steps_inner > 1.

    Splitting a step into inner steps is the same finite-difference scheme as
    using a finer grid, so (z_steps=N, z_steps_inner=2) must reproduce
    (z_steps=2N, z_steps_inner=1) at every shared z, for both Euler and AB (the
    bootstrap is the same single Euler step in both cases).
    """

    N = 10

    def _solve(self, step, z_steps, z_steps_inner):
        with open(os.path.join(JSON_DIR, "mb_solve_01.json")) as f:
            d = json.load(f)
        d["z_steps"] = z_steps
        d["z_steps_inner"] = z_steps_inner
        # Strongly absorbing, nonlinear pulse over the whole medium, so that a
        # one-inner-step offset in z is large compared with the tolerance below.
        d["z_min"], d["z_max"] = 0.0, 1.0
        d["t_min"], d["t_max"], d["t_steps"] = -0.5, 1.5, 100
        d["interaction_strengths"] = [40.0]
        d["opts"] = {"atol": 1e-11, "rtol": 1e-9}
        field = d["system"]["fields"][0]
        field["rabi_freq"] = 0.5
        field["rabi_freq_t_args"] = {"ampl": 1.0, "centre": 0.0, "fwhm": 0.2}
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "inner_steps.json")
            with open(path, "w") as f:
                json.dump(d, f)
            mbs = PropagationSolver().from_json(path)
        mbs.solve(step=step, progress=False)
        return mbs

    def test_inner_steps_match_finer_grid(self):
        for step in ("euler", "ab"):
            with self.subTest(step=step):
                coarse = self._solve(step, self.N, 2)
                fine = self._solve(step, 2 * self.N, 1)
                np.testing.assert_allclose(coarse.zlist, fine.zlist[::2])
                np.testing.assert_allclose(
                    coarse._Omegas_zt, fine._Omegas_zt[:, ::2, :], rtol=0, atol=1e-10
                )
                np.testing.assert_allclose(
                    coarse.states_zt, fine.states_zt[::2], rtol=0, atol=1e-10
                )


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


def _j_to_j_system(Jg, Je, pol, rabi_freq, rate=1.0):
    """Jg -> Je manifold with unpolarised ground state, as a system dict.

    pol "x": linear, q=+1 and q=-1 with factors -CG/sqrt(2) and +CG/sqrt(2).
    pol "sigma+": circular, q=+1 with factors CG. Decays have branching CG**2.
    """
    from math import sqrt

    from clerq.angmom import wigner_3j

    def cg(mg, me, q):
        if mg + q != me:
            return 0.0
        return (
            (-1) ** (Jg - 1 + me) * sqrt(2 * Je + 1) * wigner_3j(Jg, 1, Je, mg, q, -me)
        )

    mgs = [-Jg + i for i in range(int(round(2 * Jg)) + 1)]
    mes = [-Je + i for i in range(int(round(2 * Je)) + 1)]
    n_g = len(mgs)
    g_idx = {m: i for i, m in enumerate(mgs)}
    e_idx = {m: n_g + i for i, m in enumerate(mes)}
    qs = ((1, -1 / sqrt(2)), (-1, 1 / sqrt(2))) if pol == "x" else ((1, 1.0),)

    coupled, factors = [], []
    for mg in mgs:
        for q, s in qs:
            me = mg + q
            if me in e_idx and abs(cg(mg, me, q)) > 1e-12:
                coupled.append([g_idx[mg], e_idx[me]])
                factors.append(s * cg(mg, me, q))
    channels, decay_factors = [], []
    for me in mes:
        for mg in mgs:
            if abs(me - mg) <= 1 and cg(mg, me, me - mg) ** 2 > 1e-12:
                channels.append([g_idx[mg], e_idx[me]])
                decay_factors.append(abs(cg(mg, me, me - mg)))

    system = {
        "num_states": n_g + len(mes),
        "decays": [{"rate": rate, "channels": channels, "factors": decay_factors}],
        "initial_state": [1.0 / n_g] * n_g + [0.0] * len(mes),
        "fields": [
            {
                "coupled_levels": coupled,
                "factors": factors,
                "detuning": 0.0,
                "rabi_freq": rabi_freq,
                "rabi_freq_t_func": "square",
                "rabi_freq_t_args": {"ampl": 1.0, "on": 0.0, "off": 1.0e9},
            }
        ],
    }
    return system, np.array(factors), n_g


class TestMultiTransitionFieldSource(unittest.TestCase):
    """The source term weights each pair's coherence by its factor, like the
    Hamiltonian, so a weak resonant cw probe on a multi-transition field is
    absorbed as

        Omega(L) / Omega(0) = exp(-2 N g L sum_k(factor_k**2) / (n_g * rate))

    with n_g unpolarised ground sublevels (each excited level decays at
    2 pi rate), i.e. in proportion to the line strength sum(factor**2). It is
    not sum(factor**3), which vanishes for linear polarisation on J=1 -> 2.
    """

    N, G, L = 1.0, 1.0, 0.1

    def _transmission(self, Jg, Je, pol, step):
        system, factors, n_g = _j_to_j_system(Jg, Je, pol, rabi_freq=1.0e-3)
        mbs = PropagationSolver(
            system=system,
            t_min=0.0,
            t_max=3.0,  # many excited-state lifetimes: steady state
            t_steps=60,
            z_min=0.0,
            z_max=self.L,
            z_steps=40,
            z_steps_inner=1,
            num_density_z_func="square",
            num_density_z_args={"on": 0.0, "off": self.L, "ampl": self.N},
            interaction_strengths=[self.G],
        )
        mbs.solve(step=step)
        Omegas = mbs._Omegas_zt[0]
        ratio = (Omegas[-1, -1] / Omegas[0, -1]).real
        expected = np.exp(-2 * self.N * self.G * self.L * np.sum(factors**2) / n_g)
        return ratio, expected, factors

    def test_linear_j1_to_j2_absorbs(self):
        """J=1 -> 2 linear: sum(factor**3) = 0, but the line strength is not."""
        for step in ("euler", "ab"):
            ratio, expected, factors = self._transmission(1.0, 2.0, "x", step)
            self.assertAlmostEqual(np.sum(factors**3), 0.0, places=12)
            self.assertLess(ratio, 0.99)
            self.assertAlmostEqual(ratio, expected, delta=2e-3)

    def test_sigma_plus_j1_to_j2(self):
        """J=1 -> 2 sigma+: all factors positive, sum(f**2) != sum(f**3)."""
        for step in ("euler", "ab"):
            ratio, expected, factors = self._transmission(1.0, 2.0, "sigma+", step)
            wrong = np.exp(-2 * self.N * self.G * self.L * np.sum(factors**3) / 3)
            self.assertGreater(abs(wrong - expected), 1e-2)
            self.assertAlmostEqual(ratio, expected, delta=2e-3)

    def test_linear_j12_to_j12(self):
        ratio, expected, _ = self._transmission(0.5, 0.5, "x", "ab")
        self.assertAlmostEqual(ratio, expected, delta=2e-3)

    def test_energy_balance(self):
        """Without decay, the field energy lost per unit length equals the
        atoms' excited-state gain (a nonlinear pulse, many transitions):

            d/dz int |Omega|^2 dt = -4 N g_ang dP_e,   g_ang = 2 pi g

        This needs the source term to weight each coherence by the same
        factor as the Hamiltonian. With factor**2 it fails by 25% to 80%
        (or absorbs nothing at all for linear J=1 -> 2).
        """
        N, G, L = 1.0, 1.0, 0.3
        for Jg, Je, pol in [(1.0, 2.0, "sigma+"), (1.0, 2.0, "x"), (0.5, 0.5, "x")]:
            with self.subTest(Jg=Jg, Je=Je, pol=pol):
                system, _, n_g = _j_to_j_system(Jg, Je, pol, rabi_freq=0.4, rate=0.0)
                system["fields"][0]["rabi_freq_t_func"] = "gaussian"
                system["fields"][0]["rabi_freq_t_args"] = {
                    "ampl": 1.0,
                    "centre": 1.0,
                    "fwhm": 0.5,
                }
                mbs = PropagationSolver(
                    system=system,
                    t_min=0.0,
                    t_max=2.0,
                    t_steps=200,
                    z_min=0.0,
                    z_max=L,
                    z_steps=120,
                    z_steps_inner=1,
                    num_density_z_func="square",
                    num_density_z_args={"on": 0.0, "off": L, "ampl": N},
                    interaction_strengths=[G],
                )
                mbs.solve(step="ab")

                energy_zt = np.trapezoid(
                    np.abs(mbs._Omegas_zt[0]) ** 2, mbs.tlist, axis=1
                )
                lost = energy_zt[0] - energy_zt[-1]
                P_e = np.real(
                    np.trace(mbs.states_zt[:, :, n_g:, n_g:], axis1=2, axis2=3)
                )
                gained = (
                    4 * N * G * 2 * np.pi
                    * np.trapezoid(P_e[:, -1] - P_e[:, 0], mbs.zlist)
                )  # fmt: skip
                self.assertGreater(lost, 0.1)
                # The z-discretisation error is first order: about 1.5% here.
                self.assertAlmostEqual(gained / lost, 1.0, delta=0.03)


class TestSourceWeightsSavefile(unittest.TestCase):
    """Results solved with the old factor**2 source weights are stale unless
    all the field factors are 1, so they are not loaded."""

    def _solver(self, system, savefile):
        mbs = PropagationSolver(
            system=system,
            t_min=0.0,
            t_max=1.0,
            t_steps=20,
            z_min=0.0,
            z_max=0.1,
            z_steps=4,
            z_steps_inner=1,
            num_density_z_func="square",
            num_density_z_args={"on": 0.0, "off": 0.1, "ampl": 1.0},
            interaction_strengths=[1.0],
            savefile=savefile,
        )
        return mbs

    def _write_old_savefile(self, mbs, savefile, legacy):
        import qutip as qu

        mbs.solve()
        if legacy:
            qu.qsave((mbs._Omegas_zt, mbs.states_zt), savefile)
        else:
            meta = mbs._savefile_metadata()
            del meta["source_weights"]  # as written before the weights changed
            qu.qsave((mbs._Omegas_zt, mbs.states_zt, meta), savefile)

    def test_current_savefile_loads(self):
        import tempfile

        system, _, _ = _j_to_j_system(0.5, 0.5, "x", rabi_freq=1.0e-3)
        with tempfile.TemporaryDirectory() as tmpdir:
            mbs = self._solver(system, os.path.join(tmpdir, "result"))
            mbs.solve()
            mbs.load_results()

    def test_old_savefile_with_nonunit_factors_rejected(self):
        import tempfile

        system, factors, _ = _j_to_j_system(0.5, 0.5, "x", rabi_freq=1.0e-3)
        self.assertTrue(np.any(factors != 1.0))
        for legacy in (False, True):
            with self.subTest(legacy=legacy), tempfile.TemporaryDirectory() as tmpdir:
                savefile = os.path.join(tmpdir, "result")
                mbs = self._solver(system, savefile)
                self._write_old_savefile(mbs, savefile, legacy)
                with self.assertRaises(ValueError):
                    mbs.load_results()
                # solve() recalculates rather than using the stale results,
                # and writes a savefile that loads.
                mbs.solve(recalc=False)
                mbs.load_results()

    def test_old_savefile_with_unit_factors_loads(self):
        import tempfile

        system = {
            "num_states": 2,
            "decays": [{"rate": 1.0, "channels": [[0, 1]], "factors": [1.0]}],
            "initial_state": [1.0, 0.0],
            "fields": [
                {
                    "coupled_levels": [[0, 1]],
                    "factors": [1.0],
                    "detuning": 0.0,
                    "rabi_freq": 1.0e-3,
                }
            ],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            savefile = os.path.join(tmpdir, "result")
            mbs = self._solver(system, savefile)
            self._write_old_savefile(mbs, savefile, legacy=False)
            mbs.load_results()


class TestRabiFreqsZt(unittest.TestCase):
    """rabi_freqs_zt is the public cyclic view of the internal angular field."""

    def _solver(self):
        mbs = PropagationSolver().from_json(os.path.join(JSON_DIR, "mb_solve_01.json"))
        mbs.solve(recalc=True, progress=False)
        return mbs

    def test_cyclic_view_of_angular_field(self):
        mbs = self._solver()
        np.testing.assert_allclose(mbs.rabi_freqs_zt, mbs._Omegas_zt / (2 * np.pi))

    def test_input_slice_is_the_problem_rabi_freq(self):
        mbs = self._solver()
        f = mbs.system.fields[0]
        expected = f.rabi_freq * f.rabi_freq_t_func(mbs.tlist, f.rabi_freq_t_args)
        np.testing.assert_allclose(mbs.rabi_freqs_zt[0, 0], expected)


class TestDriveConvention(unittest.TestCase):
    """The savefile and JSON record which drive (complex or Re(Omega)) was used."""

    def _solver(self, savefile=None, detuning=0.0, real_drive=False, classes=None):
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            return PropagationSolver(
                system={
                    "num_states": 2,
                    "decays": [{"channels": [[0, 1]], "rate": 1.0}],
                    "fields": [
                        {
                            "coupled_levels": [[0, 1]],
                            "detuning": detuning,
                            "rabi_freq": 1.0e-3,
                            "rabi_freq_t_func": "gaussian",
                            "rabi_freq_t_args": {
                                "ampl": 1.0,
                                "centre": 0.5,
                                "fwhm": 0.3,
                            },
                        }
                    ],
                },
                t_min=0.0,
                t_max=1.0,
                t_steps=20,
                z_max=0.1,
                z_steps=3,
                z_steps_inner=1,
                interaction_strengths=[1.0],
                velocity_classes=classes,
                savefile=savefile,
                real_drive=real_drive,
            )

    def _write(self, mbs, savefile, drive):
        import qutip as qu

        mbs.solve(check_stability=False)
        meta = mbs._savefile_metadata()
        if drive is None:
            del meta["drive"]  # as written before the marker existed
        else:
            meta["drive"] = drive
        qu.qsave((mbs._Omegas_zt, mbs.states_zt, meta), savefile)

    def test_json_round_trip_keeps_real_drive(self):
        mbs = self._solver(real_drive=True)
        self.assertTrue(mbs.get_json_dict()["real_drive"])
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            again = PropagationSolver.from_json_str(mbs.to_json_str())
        self.assertTrue(again.real_drive)
        self.assertEqual(again.to_json_str(), mbs.to_json_str())

    def test_default_json_has_no_real_drive_key(self):
        mbs = self._solver()
        self.assertNotIn("real_drive", mbs.get_json_dict())
        self.assertFalse(PropagationSolver.from_json_str(mbs.to_json_str()).real_drive)

    def test_current_savefile_loads(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            f = os.path.join(tmpdir, "r")
            mbs = self._solver(f, detuning=0.5)
            mbs.solve(check_stability=False)
            mbs.load_results()
            self.assertEqual(mbs._savefile_metadata()["drive"], "complex")

    def test_old_savefile_rejected_when_detuned(self):
        import tempfile

        for drive in (None, "real"):
            with tempfile.TemporaryDirectory() as tmpdir:
                f = os.path.join(tmpdir, "r")
                mbs = self._solver(f, detuning=0.5)
                self._write(mbs, f, drive)
                with self.assertRaisesRegex(ValueError, "drive"):
                    mbs.load_results()

    def test_old_savefile_rejected_for_asymmetric_classes(self):
        import tempfile

        classes = {
            "thermal_delta_min": 0.0,
            "thermal_delta_max": 2.0,
            "thermal_delta_steps": 2,
            "thermal_width": 1.0,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            f = os.path.join(tmpdir, "r")
            mbs = self._solver(f, classes=classes)
            self._write(mbs, f, None)
            with self.assertRaisesRegex(ValueError, "drive"):
                mbs.load_results()

    def test_old_resonant_symmetric_savefile_still_loads(self):
        import tempfile

        classes = {
            "thermal_delta_min": -2.0,
            "thermal_delta_max": 2.0,
            "thermal_delta_steps": 4,
            "thermal_width": 1.0,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            f = os.path.join(tmpdir, "r")
            mbs = self._solver(f, classes=classes)
            self._write(mbs, f, None)
            mbs.load_results()

    def test_real_drive_savefile_rejected_by_complex_solver(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            f = os.path.join(tmpdir, "r")
            old = self._solver(f, detuning=0.5, real_drive=True)
            old.solve(check_stability=False)
            self.assertEqual(old._savefile_metadata()["drive"], "real")
            # same savefile, default drive: the hash (real_drive key) differs
            with self.assertRaises(ValueError):
                self._solver(f, detuning=0.5).load_results()
