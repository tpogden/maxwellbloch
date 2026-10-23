"""
Unit tests for the DC B-field (Zeeman shift) support.

The Zeeman shift of a sublevel is E = +g_F m_F B in natural units (B stands for
mu_B B / h in the frequency unit of the system's energies).
"""

import copy
import json

import numpy as np
import pytest
from scipy.signal import find_peaks

from clerq import PropagationSolver, System, hyperfine, spectral


def two_level_problem(b_field=None, fwhm=0.05):
    """Weak broadband probe on a two-level system, spectrally resolvable."""
    d = {
        "system": {
            "num_states": 2,
            "decays": [{"channels": [[0, 1]], "rate": 1.0}],
            "fields": [
                {
                    "coupled_levels": [[0, 1]],
                    "rabi_freq": 1e-3,
                    "rabi_freq_t_func": "gaussian",
                    "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": fwhm},
                }
            ],
        },
        "t_min": -2.0,
        "t_max": 30.0,
        "t_steps": 3200,
        "z_min": -0.2,
        "z_max": 1.2,
        "z_steps": 20,
        "z_steps_inner": 1,
        "interaction_strengths": [0.5],
    }
    if b_field is not None:
        d["system"]["b_field"] = b_field
    return d


def solve(problem):
    solver = PropagationSolver.from_json_str(json.dumps(problem))
    solver.solve(progress=False)
    return solver


def absorption_peaks(solver, min_height_frac=0.3):
    freqs = spectral.freq_list(solver)
    absorb = spectral.absorption(solver, 0, -1)
    idx, _ = find_peaks(absorb, height=min_height_frac * absorb.max())
    return freqs, absorb, freqs[idx]


class TestGF:
    def test_lande_g_J_rb87(self):
        assert hyperfine.lande_g_J(L=0, J=0.5) == pytest.approx(2.0023193)
        assert hyperfine.lande_g_J(L=1, J=1.5) == pytest.approx(1.3362, abs=3e-3)
        assert hyperfine.lande_g_J(L=1, J=0.5) == pytest.approx(0.6659, abs=3e-3)

    def test_g_F_rb87_ground(self):
        """Steck: g_F = +1/2 for F=2 and -1/2 for F=1 of 5S1/2."""
        g_I = -0.0009951
        f2 = hyperfine.LevelF(I=1.5, J=0.5, F=2, L=0, g_I=g_I)
        f1 = hyperfine.LevelF(I=1.5, J=0.5, F=1, L=0, g_I=g_I)
        assert f2.g_F == pytest.approx(0.4998, abs=1e-4)
        assert f1.g_F == pytest.approx(-0.5018, abs=1e-4)

    def test_g_F_excited(self):
        p32 = hyperfine.LevelF(I=1.5, J=1.5, F=3, L=1, g_I=0.0)
        assert p32.g_F == pytest.approx(2.0 / 3.0, abs=1e-3)
        p12 = hyperfine.LevelF(I=1.5, J=0.5, F=1, L=1, g_I=0.0)
        assert p12.g_F == pytest.approx(-1.0 / 6.0, abs=1e-3)

    def test_explicit_g_F_wins_and_missing_is_none(self):
        assert hyperfine.LevelF(I=1.5, J=0.5, F=1, g_F=0.25, L=0).g_F == 0.25
        assert hyperfine.LevelF(I=1.5, J=0.5, F=1).g_F is None

    def test_missing_g_I_warns_and_neglects_nuclear_term(self):
        with pytest.warns(UserWarning, match="g_I not given"):
            f1 = hyperfine.LevelF(I=1.5, J=0.5, F=1, L=0)
        assert f1.g_F == pytest.approx(-0.5006, abs=1e-4)

    def test_f0_has_zero_g_F(self):
        assert hyperfine.LevelF(I=1.5, J=1.5, F=0, L=1, g_I=0.0).g_F == 0.0

    def test_atom_zeeman_levels(self):
        atom = hyperfine.Atom1e(element="Rb", isotope="87")
        atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1, g_F=-0.5))
        atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1, g_F=1 / 6))
        z = atom.get_zeeman_levels()
        assert z["m_F"] == [-1.0, 0.0, 1.0] * 2
        assert z["g_F"] == [-0.5] * 3 + [1 / 6] * 3

    def test_atom_missing_g_F_raises(self):
        atom = hyperfine.Atom1e()
        atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))
        with pytest.raises(ValueError, match="no g_F"):
            atom.get_zeeman_levels()

    def test_json_dict_includes_g_F_only_when_set(self):
        assert "g_F" not in hyperfine.LevelF(I=1.5, J=0.5, F=1).get_json_dict()
        assert (
            hyperfine.LevelF(I=1.5, J=0.5, F=1, g_F=0.5).get_json_dict()["g_F"] == 0.5
        )

    def test_gauss_to_natural(self):
        assert hyperfine.gauss_to_natural(1.0) == pytest.approx(1.399624493)


class TestSystemBField:
    BF = {"magnitude": 2.0, "g_F": [0.5, -1.0, 2.0], "m_F": [1.0, 1.0, -0.5]}

    def test_H_0_diagonal(self):
        s = System(num_states=3, energies=[0.0, 1.0, 5.0], b_field=self.BF)
        expected = (
            2 * np.pi * (np.array([0.0, 1.0, 5.0]) + 2.0 * np.array([0.5, -1, -1]))
        )
        np.testing.assert_allclose(np.diag(s.H_0.full()).real, expected)
        np.testing.assert_allclose(s.H_0.full() - np.diag(np.diag(s.H_0.full())), 0)

    def test_gauss_units(self):
        bf = {**self.BF, "magnitude": 3.0, "units": "gauss"}
        s = System(num_states=3, b_field=bf)
        nat = System(num_states=3, b_field={**self.BF, "magnitude": 3 * 1.399624493})
        np.testing.assert_allclose(s.H_0.full(), nat.H_0.full())

    def test_zero_field_identical(self):
        plain = System(num_states=3, energies=[0.0, 1.0, 5.0])
        zero = System(
            num_states=3, energies=[0.0, 1.0, 5.0], b_field={**self.BF, "magnitude": 0}
        )
        assert np.array_equal(plain.H_0.full(), zero.H_0.full())

    def test_json_round_trip(self):
        s = System(num_states=3, b_field=self.BF)
        s2 = System.from_json_str(s.to_json_str())
        assert s2.b_field == s.b_field
        assert "b_field" not in System(num_states=3).get_json_dict()

    @pytest.mark.parametrize(
        "bf, msg",
        [
            ({"g_F": [1, 1, 1]}, "'m_F'"),
            ({"g_F": [1, 1], "m_F": [1, 1, 1]}, "'g_F'"),
            ({**BF, "units": "tesla"}, "units"),
            ({**BF, "direction": [0, 0, 1]}, "Unknown"),
        ],
    )
    def test_validation(self, bf, msg):
        with pytest.raises(ValueError, match=msg):
            System(num_states=3, b_field=bf)


class TestZeroFieldRegression:
    def test_propagation_unchanged_at_B0(self):
        base = two_level_problem()
        with_b = two_level_problem(
            {"magnitude": 0.0, "g_F": [1.0, 1.5], "m_F": [1.0, 1.0]}
        )
        a, b = solve(base), solve(with_b)
        assert np.array_equal(a._Omegas_zt, b._Omegas_zt)

    def test_hash_unchanged_without_b_field(self):
        """No block: the JSON (hence the .qu hash) has no b_field and round-trips;
        adding a block changes it."""
        plain = PropagationSolver.from_json_str(json.dumps(two_level_problem()))
        assert "b_field" not in plain.to_json_str()
        again = PropagationSolver.from_json_str(plain.to_json_str())
        assert again.to_json_str() == plain.to_json_str()
        with_b = PropagationSolver.from_json_str(
            json.dumps(
                two_level_problem({"magnitude": 0.0, "g_F": [1, 1], "m_F": [1, 1]})
            )
        )
        assert with_b.to_json_str() != plain.to_json_str()

    def test_editing_b_field_rebuilds_H_0(self):
        s = System(num_states=2, b_field={"g_F": [1, 1], "m_F": [1, 2]})
        s.b_field["magnitude"] = 9.0  # a copy: no effect
        np.testing.assert_allclose(s.H_0.full(), 0)
        s.b_field = {"magnitude": 1.0, "g_F": [1, 1], "m_F": [1, 2]}
        np.testing.assert_allclose(
            np.diag(s.H_0.full()).real, 2 * np.pi * np.array([1, 2])
        )
        assert s.get_json_dict()["b_field"]["magnitude"] == 1.0


class TestSpectrum:
    @pytest.mark.parametrize(
        "B, ge, gg", [(1.0, 1.5, -0.5), (-1.0, 1.5, -0.5), (2.0, 0.5, 1.0)]
    )
    def test_two_level_line_shift(self, B, ge, gg):
        """The line moves by (g_e m_e - g_g m_g) B."""
        bf = {"magnitude": B, "g_F": [gg, ge], "m_F": [1.0, 1.0]}
        _, absorb, peaks = absorption_peaks(solve(two_level_problem(bf)))
        assert len(peaks) == 1
        assert peaks[0] == pytest.approx((ge - gg) * B, abs=0.05)

    def test_gauss_matches_natural(self):
        nat = {"magnitude": 1.399624493, "g_F": [0.0, 1.0], "m_F": [0.0, 1.0]}
        gauss = {**nat, "magnitude": 1.0, "units": "gauss"}
        a = solve(two_level_problem(nat))._Omegas_zt
        b = solve(two_level_problem(gauss))._Omegas_zt
        np.testing.assert_allclose(a, b, rtol=1e-10)


def hyperfine_problem(q, B, fwhm=0.05):
    """Rb87 5S1/2 F=1 -> 5P1/2 F'=1 with a weak probe of polarisation q."""
    atom = hyperfine.Atom1e(element="Rb", isotope="87")
    g_I = -0.0009951
    atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1, L=0, g_I=g_I))
    atom.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1, L=1, g_I=g_I))
    problem = two_level_problem(fwhm=fwhm)
    coupled = atom.get_coupled_levels([0], [1])
    problem["system"].update(
        {
            "num_states": 6,
            "decays": [
                {
                    "channels": coupled,
                    "rate": 1.0,
                    "factors": atom.get_decay_factors([0], [1]).tolist(),
                }
            ],
            "initial_state": [1 / 3] * 3 + [0.0] * 3,
            "b_field": {"magnitude": B, **atom.get_zeeman_levels()},
        }
    )
    problem["system"]["fields"][0].update(
        {
            "coupled_levels": coupled,
            "factors": atom.get_clebsch_hf_factors([0], [1], q=q).tolist(),
        }
    )
    return problem, atom


def expected_lines(atom, q, B):
    """Line positions (g' m' - g m) B of the coupled sublevel pairs with q."""
    mF = [lvl["mF"] for lvl in atom.get_mF_list()]
    g_F = atom.get_g_F_list()
    pairs = atom.get_coupled_levels([0], [1])
    factors = atom.get_clebsch_hf_factors([0], [1], q=q)
    return sorted(
        (g_F[b] * mF[b] - g_F[a] * mF[a]) * B
        for (a, b), f in zip(pairs, factors, strict=True)
        if abs(f) > 1e-12
    )


def assert_lines(solver, expected, window=1.0):
    """The spectrum has a peak at each expected line and no other absorption.

    The natural-linewidth Lorentzian wings of the lines are physical, so "no
    other absorption" starts ``window`` away from each line.
    """
    freqs = spectral.freq_list(solver)
    absorb = spectral.absorption(solver, 0, -1)
    inband = np.abs(freqs) < 20
    top = absorb[inband].max()
    covered = np.zeros_like(inband)
    for line in expected:
        near = np.abs(freqs - line) < window
        covered |= near
        assert absorb[near].max() > 0.05 * top
        assert abs(freqs[near][np.argmax(absorb[near])] - line) < 0.5
    assert absorb[inband & ~covered].max() < 0.02 * top


class TestPolarisation:
    B = 20.0

    @pytest.mark.parametrize("q", [1, -1])
    def test_sublevel_resolved_lines(self, q):
        """Each allowed sublevel pair has its line at (g' m' - g m) B."""
        problem, atom = hyperfine_problem(q=q, B=self.B)
        expected = expected_lines(atom, q, self.B)
        assert len(expected) == 2
        assert_lines(solve(problem), expected)

    def test_sigma_plus_and_minus_mirror(self):
        """Flipping q and the sign of B gives the same spectrum."""
        p_plus, _ = hyperfine_problem(q=1, B=self.B)
        p_minus, _ = hyperfine_problem(q=-1, B=-self.B)
        a = spectral.absorption(solve(p_plus), 0, -1)
        b = spectral.absorption(solve(p_minus), 0, -1)
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-8 * np.abs(a).max())

    def test_sigma_plus_and_minus_lines_opposite(self):
        """At the same B the sigma- lines are the sigma+ lines reflected."""
        _, atom = hyperfine_problem(q=1, B=self.B)
        plus = expected_lines(atom, 1, self.B)
        minus = expected_lines(atom, -1, self.B)
        np.testing.assert_allclose(sorted(minus), sorted(-np.array(plus)), atol=1e-9)


class TestAutodiff:
    @pytest.fixture(autouse=True)
    def _jax(self):
        jax = pytest.importorskip("jax")
        pytest.importorskip("diffrax")
        jax.config.update("jax_enable_x64", True)

    def problem(self, B):
        d = copy.deepcopy(
            two_level_problem(
                {"magnitude": B, "g_F": [-0.5, 1.5], "m_F": [1.0, 1.0]}, fwhm=0.4
            )
        )
        d["t_min"], d["t_max"], d["t_steps"] = -1.0, 4.0, 60
        d["z_steps"] = 4
        d["system"]["fields"][0]["rabi_freq"] = 0.5
        d["interaction_strengths"] = [2.0]
        return d

    def solver(self, B):
        s = PropagationSolver.from_json_str(json.dumps(self.problem(B)))
        s.build_opts({"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 10**6})
        return s

    def test_matches_solve_with_nonzero_B(self):
        s = self.solver(0.7)
        s.solve(progress=False)
        got = np.asarray(s.to_autodiff().propagate())
        np.testing.assert_allclose(got, s.rabi_freqs_zt, rtol=1e-4, atol=1e-7)

    def test_b_field_param_matches_rebuilt_system(self):
        """Passing b_field as a parameter equals building the system at that B."""
        ad = self.solver(0.0).to_autodiff()
        s = self.solver(0.7)
        ref = np.asarray(s.to_autodiff().propagate())
        got = np.asarray(ad.propagate({"b_field": 0.7}))
        np.testing.assert_allclose(got, ref, rtol=1e-8, atol=1e-10)

    def test_default_b_field_param_is_systems(self):
        ad = self.solver(0.7).to_autodiff()
        assert ad.default_params()["b_field"] == pytest.approx(0.7)

    def test_gradient_wrt_B_matches_finite_difference(self):
        import jax.numpy as jnp

        ad = self.solver(0.7).to_autodiff()

        def loss(omega):
            return jnp.sum(jnp.abs(omega[0, -1]) ** 2)

        _, grads = ad.value_and_grad(loss, {"b_field": 0.7})
        h = 1e-4
        up = float(loss(ad.propagate({"b_field": 0.7 + h})))
        dn = float(loss(ad.propagate({"b_field": 0.7 - h})))
        assert float(grads["b_field"]) == pytest.approx((up - dn) / (2 * h), rel=1e-4)

    def test_b_field_param_without_block_raises(self):
        d = self.problem(0.7)
        del d["system"]["b_field"]
        ad = PropagationSolver.from_json_str(json.dumps(d)).to_autodiff()
        assert "b_field" not in ad.default_params()
        with pytest.raises(ValueError, match="no b_field block"):
            ad.propagate({"b_field": 0.7})

    def test_zero_zeeman_block_keeps_b_field_param(self):
        d = self.problem(0.7)
        d["system"]["b_field"]["m_F"] = [0.0, 0.0]
        ad = PropagationSolver.from_json_str(json.dumps(d)).to_autodiff()
        assert ad.default_params()["b_field"] == pytest.approx(0.7)

    def test_gauss_units_default_param_is_natural(self):
        d = self.problem(0.7)
        d["system"]["b_field"].update(magnitude=0.5, units="gauss")
        ad = PropagationSolver.from_json_str(json.dumps(d)).to_autodiff()
        nat = 0.5 * hyperfine.MU_B_OVER_H_MHZ_PER_GAUSS
        assert ad.default_params()["b_field"] == pytest.approx(nat)
        ref = np.asarray(self.solver(nat).to_autodiff().propagate())
        np.testing.assert_allclose(
            np.asarray(ad.propagate()), ref, rtol=1e-8, atol=1e-10
        )
