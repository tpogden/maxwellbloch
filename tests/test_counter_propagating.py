"""Tests for counter-propagating field support (A0).

Physics cases:
1. Undepleted regime: depletion check passes silently.
2. Moderately depleted regime: UserWarning emitted.
3. Heavily depleted regime: CounterPropagatingDepletionError raised.
4. Doppler sign correctness: counter-propagating coupling gives Doppler-free
   two-photon resonance (narrower EIT linewidth than co-propagating).
5. Regression: omitting new attrs gives bit-exact output vs existing fixtures.

Thomas Ogden <t@ogden.eu>
"""

import unittest
import warnings

import numpy as np

from maxwellbloch import mb_solve
from maxwellbloch.exceptions import CounterPropagatingDepletionError
from maxwellbloch.field import Field

# ---------------------------------------------------------------------------
# Field schema / init tests
# ---------------------------------------------------------------------------


class TestFieldCounterPropInit(unittest.TestCase):
    def test_defaults(self):
        """Default field has counter_propagating=False, factor_doppler_shift=1.0."""
        f = Field(coupled_levels=[[0, 1]])
        self.assertFalse(f.counter_propagating)
        self.assertEqual(f.factor_doppler_shift, 1.0)

    def test_counter_propagating_sets_factor(self):
        """counter_propagating=True sets factor_doppler_shift to -1.0."""
        f = Field(coupled_levels=[[0, 1]], counter_propagating=True)
        self.assertTrue(f.counter_propagating)
        self.assertEqual(f.factor_doppler_shift, -1.0)

    def test_factor_doppler_shift_negative(self):
        """factor_doppler_shift can be set to a negative float directly."""
        f = Field(coupled_levels=[[0, 1]], factor_doppler_shift=-0.62)
        self.assertFalse(f.counter_propagating)
        self.assertAlmostEqual(f.factor_doppler_shift, -0.62)

    def test_counter_propagating_false_with_negative_factor(self):
        """counter_propagating=False and factor_doppler_shift=-1.0 is legal."""
        f = Field(
            coupled_levels=[[0, 1]],
            counter_propagating=False,
            factor_doppler_shift=-1.0,
        )
        self.assertFalse(f.counter_propagating)
        self.assertEqual(f.factor_doppler_shift, -1.0)

    def test_conflicting_flags_raise(self):
        """counter_propagating=True with factor_doppler_shift>0 raises ValueError."""
        with self.assertRaises(ValueError):
            Field(
                coupled_levels=[[0, 1]],
                counter_propagating=True,
                factor_doppler_shift=1.0,
            )

    def test_json_roundtrip(self):
        """New attrs survive JSON serialisation round-trip."""
        f = Field(
            coupled_levels=[[0, 1]],
            counter_propagating=True,
            label="coupling",
        )
        f2 = Field.from_json_str(f.to_json_str())
        self.assertTrue(f2.counter_propagating)
        self.assertEqual(f2.factor_doppler_shift, -1.0)

    def test_json_roundtrip_factor_only(self):
        """factor_doppler_shift=-0.62 survives JSON round-trip."""
        f = Field(coupled_levels=[[0, 1]], factor_doppler_shift=-0.62)
        f2 = Field.from_json_str(f.to_json_str())
        self.assertAlmostEqual(f2.factor_doppler_shift, -0.62)


# ---------------------------------------------------------------------------
# Regression: bit-exact output when new attrs are absent from JSON
# ---------------------------------------------------------------------------

# Minimal 3-level EIT config reused across depletion tests. Parameters chosen
# so the coupling field is strong and barely depletes (OD ≈ 0.5 per unit length).
_LADDER_BASE_JSON = """
{
  "atom": {
    "num_states": 3,
    "decays": [{"channels": [[0, 1], [1, 2]], "rate": 1.0}],
    "fields": [
      {
        "label": "probe",
        "coupled_levels": [[0, 1]],
        "rabi_freq": 1.0,
        "rabi_freq_t_func": "sech",
        "rabi_freq_t_args": {"ampl": 0.5, "centre": 0.0, "width": 0.5}
      },
      {
        "label": "coupling",
        "coupled_levels": [[1, 2]],
        "rabi_freq": 5.0,
        "rabi_freq_t_func": "ramp_onoff",
        "rabi_freq_t_args": {"ampl": 1.0, "fwhm": 0.2, "on": -2.0, "off": 8.0},
        "counter_propagating": true
      }
    ]
  },
  "t_min": -2.0,
  "t_max": 8.0,
  "t_steps": 50,
  "z_min": 0.0,
  "z_max": 1.0,
  "z_steps": 5,
  "z_steps_inner": 2,
  "interaction_strengths": [0.5, 0.5]
}
"""

# Same config but without counter_propagating on the coupling. Used as the
# regression baseline — output must be identical when factor_doppler_shift=1
# (no velocity classes, so Doppler shift is irrelevant).
_LADDER_NO_CP_JSON = """
{
  "atom": {
    "num_states": 3,
    "decays": [{"channels": [[0, 1], [1, 2]], "rate": 1.0}],
    "fields": [
      {
        "label": "probe",
        "coupled_levels": [[0, 1]],
        "rabi_freq": 1.0,
        "rabi_freq_t_func": "sech",
        "rabi_freq_t_args": {"ampl": 0.5, "centre": 0.0, "width": 0.5}
      },
      {
        "label": "coupling",
        "coupled_levels": [[1, 2]],
        "rabi_freq": 5.0,
        "rabi_freq_t_func": "ramp_onoff",
        "rabi_freq_t_args": {"ampl": 1.0, "fwhm": 0.2, "on": -2.0, "off": 8.0}
      }
    ]
  },
  "t_min": -2.0,
  "t_max": 8.0,
  "t_steps": 50,
  "z_min": 0.0,
  "z_max": 1.0,
  "z_steps": 5,
  "z_steps_inner": 2,
  "interaction_strengths": [0.5, 0.5]
}
"""


class TestRegression(unittest.TestCase):
    def test_no_velocity_classes_bit_exact(self):
        """Without velocity classes, counter_propagating has no effect on Omegas_zt.

        With no Doppler broadening (no velocity_classes), the thermal_delta_list
        contains only [0.0], so factor_doppler_shift multiplies zero — output is
        bit-exact regardless of the flag.
        """
        mbs_cp = mb_solve.MBSolve.from_json_str(_LADDER_BASE_JSON)
        mbs_no_cp = mb_solve.MBSolve.from_json_str(_LADDER_NO_CP_JSON)

        mbs_cp.mbsolve(recalc=True, check_counter_prop_depletion=False)
        mbs_no_cp.mbsolve(recalc=True, check_counter_prop_depletion=False)

        np.testing.assert_array_equal(mbs_cp.Omegas_zt, mbs_no_cp.Omegas_zt)
        np.testing.assert_array_equal(mbs_cp.states_zt, mbs_no_cp.states_zt)


# ---------------------------------------------------------------------------
# Depletion check tests
# ---------------------------------------------------------------------------

# Two-level absorber with a weak Gaussian probe. At high OD the probe depletes
# substantially; at low OD it passes through. We mark the probe field as
# counter_propagating so the depletion check fires.
# Interaction-strength values that hit each depletion band come from numerical
# calibration (see plans/phase-a/A0-counter-propagating-field-support.md):
#   strength=0.003 → ~0.2%  (undepleted, < 1%)
#   strength=0.07  → ~5.0%  (moderate,   1–10%)
#   strength=0.5   → ~30.5% (heavy,       > 10%)


def _two_level_absorber_json(interaction_strength: float) -> str:
    return f"""
{{
  "atom": {{
    "num_states": 2,
    "decays": [{{"channels": [[0, 1]], "rate": 1.0}}],
    "fields": [
      {{
        "label": "probe",
        "coupled_levels": [[0, 1]],
        "rabi_freq": 0.001,
        "rabi_freq_t_func": "gaussian",
        "rabi_freq_t_args": {{"ampl": 1.0, "centre": 0.0, "fwhm": 1.0}},
        "counter_propagating": true
      }}
    ]
  }},
  "t_min": -2.0,
  "t_max": 10.0,
  "t_steps": 120,
  "z_min": 0.0,
  "z_max": 1.0,
  "z_steps": 30,
  "z_steps_inner": 2,
  "interaction_strengths": [{interaction_strength}]
}}
"""


class TestDepletionCheckUndepleted(unittest.TestCase):
    def test_no_warning_at_low_od(self):
        """Weak probe at very low OD: depletion < 1%, check passes silently."""
        mbs = mb_solve.MBSolve.from_json_str(_two_level_absorber_json(0.003))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mbs.mbsolve(recalc=True)
        depletion_warnings = [
            x
            for x in w
            if issubclass(x.category, UserWarning)
            and "depleted" in str(x.message).lower()
        ]
        self.assertEqual(len(depletion_warnings), 0)

        # Confirm depletion really is < 1%
        Omega = mbs.Omegas_zt[0]
        peak_zmin = np.max(np.abs(Omega[0]))
        peak_zmax = np.max(np.abs(Omega[-1]))
        depletion = 1.0 - peak_zmax / peak_zmin
        self.assertLess(depletion, 0.01)


class TestDepletionCheckModerate(unittest.TestCase):
    def test_warning_at_moderate_od(self):
        """Probe at moderate OD: depletion between 1% and 10%, UserWarning emitted."""
        mbs = mb_solve.MBSolve.from_json_str(_two_level_absorber_json(0.07))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mbs.mbsolve(recalc=True)
        depletion_warnings = [
            x
            for x in w
            if issubclass(x.category, UserWarning)
            and "depleted" in str(x.message).lower()
        ]
        self.assertGreater(len(depletion_warnings), 0)
        # Warning text should mention "forward-propagation"
        self.assertIn(
            "forward-propagation",
            str(depletion_warnings[0].message).lower(),
        )

        # Confirm depletion is in the 1–10% band
        Omega = mbs.Omegas_zt[0]
        peak_zmin = np.max(np.abs(Omega[0]))
        peak_zmax = np.max(np.abs(Omega[-1]))
        depletion = 1.0 - peak_zmax / peak_zmin
        self.assertGreaterEqual(depletion, 0.01)
        self.assertLess(depletion, 0.10)


class TestDepletionCheckHeavy(unittest.TestCase):
    def test_error_at_high_od(self):
        """Probe at high OD: depletion > 10%, CounterPropagatingDepletionError raised."""
        mbs = mb_solve.MBSolve.from_json_str(_two_level_absorber_json(0.5))
        with self.assertRaises(CounterPropagatingDepletionError):
            mbs.mbsolve(recalc=True)

    def test_suppressed_with_flag(self):
        """check_counter_prop_depletion=False suppresses the error."""
        mbs = mb_solve.MBSolve.from_json_str(_two_level_absorber_json(0.5))
        # Should not raise
        mbs.mbsolve(recalc=True, check_counter_prop_depletion=False)


# ---------------------------------------------------------------------------
# Doppler sign correctness
# ---------------------------------------------------------------------------

# Direct sign test using a 2-level system.
#
# A probe is detuned by D from the atomic resonance. A single velocity class
# is chosen so that its Doppler shift exactly cancels the probe detuning for
# the COUNTER-propagating case, while doubling it for the CO-propagating case.
#
# Convention in MBSolve: thermal_delta_list[j] = 2π * thermal_delta_value.
# This is added directly to f.detuning, which _build_H_Delta multiplies by 2π.
# So the resonance condition for counter-prop is:
#     f.detuning + (-1) * 2π * delta = 0  →  delta = D / (2π)
#
# At delta = D/(2π) ≈ 0.159 * D:
#   co-prop:      eff_detuning ≈ D + D = 2D  →  off resonance → low absorption
#   counter-prop: eff_detuning ≈ D - D = 0   →  on resonance  → high absorption
#
# Numerically calibrated: counter-prop gives ~11x more absorption than co-prop.

_DOPPLER_SIGN_D = 1.0  # probe detuning in γ units
_DOPPLER_SIGN_DELTA = _DOPPLER_SIGN_D / (2.0 * np.pi)  # single velocity class value


def _doppler_sign_json(counter_propagating: bool) -> str:
    cp_str = "true" if counter_propagating else "false"
    d = _DOPPLER_SIGN_DELTA
    return f"""
{{
  "atom": {{
    "num_states": 2,
    "decays": [{{"channels": [[0, 1]], "rate": 1.0}}],
    "fields": [
      {{
        "label": "probe",
        "coupled_levels": [[0, 1]],
        "detuning": {_DOPPLER_SIGN_D},
        "rabi_freq": 0.001,
        "rabi_freq_t_func": "gaussian",
        "rabi_freq_t_args": {{"ampl": 1.0, "centre": 0.0, "fwhm": 1.0}},
        "counter_propagating": {cp_str}
      }}
    ]
  }},
  "t_min": -2.0,
  "t_max": 8.0,
  "t_steps": 80,
  "z_min": 0.0,
  "z_max": 1.0,
  "z_steps": 5,
  "z_steps_inner": 2,
  "interaction_strengths": [0.1],
  "velocity_classes": {{
    "thermal_width": 10.0,
    "thermal_delta_min": {d},
    "thermal_delta_max": {d},
    "thermal_delta_steps": 0,
    "thermal_delta_inner_min": {d},
    "thermal_delta_inner_max": {d},
    "thermal_delta_inner_steps": 0
  }}
}}
"""


class TestDopplerSignCorrectness(unittest.TestCase):
    def test_counter_prop_increases_absorption_at_resonance(self):
        """Counter-propagating flag flips the Doppler sign, bringing probe on resonance.

        With probe detuning D and a single velocity class at delta = D/(2π):
        - Co-propagating:      effective detuning = D + 2π*(D/2π) =  2D (off resonance)
        - Counter-propagating: effective detuning = D - 2π*(D/2π) =  0  (on resonance)

        Counter-propagating should therefore give significantly more absorption.
        Numerically calibrated: counter-prop depletion / co-prop depletion ≈ 11×.
        """
        mbs_co = mb_solve.MBSolve.from_json_str(_doppler_sign_json(False))
        mbs_cp = mb_solve.MBSolve.from_json_str(_doppler_sign_json(True))

        mbs_co.mbsolve(recalc=True, check_counter_prop_depletion=False)
        mbs_cp.mbsolve(recalc=True, check_counter_prop_depletion=False)

        def _depletion(mbs: mb_solve.MBSolve) -> float:
            Omega = mbs.Omegas_zt[0]
            peak_zmin = np.max(np.abs(Omega[0]))
            peak_zmax = np.max(np.abs(Omega[-1]))
            return 1.0 - peak_zmax / peak_zmin

        dep_co = _depletion(mbs_co)
        dep_cp = _depletion(mbs_cp)

        # Counter-prop should give substantially more absorption (at least 3×)
        self.assertGreater(
            dep_cp,
            3.0 * dep_co,
            msg=(
                f"Counter-propagating (on resonance) depletion {dep_cp:.4%} "
                f"should be at least 3× co-propagating (off resonance) {dep_co:.4%}. "
                "This tests that the Doppler sign flips correctly."
            ),
        )


if __name__ == "__main__":
    unittest.main()
