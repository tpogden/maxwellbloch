# -*- coding: utf-8 -*-

"""Comprehensive 'benchmax' benchmark — single config, all major features.

System: Rb 87 V-simulton on D1 + D2, full hyperfine, thermal atoms.

  5P₁/₂ F'=1  (states 3-5)   ← D1 upper (795 nm)
       |
  [field 0: sech, 2π]
       |
  5S₁/₂ F=1   (states 0-2)
       |
  [field 1: sech, 2π]
       |
  5P₃/₂ F'=2  (states 6-10)  ← D2 upper (780 nm)

Both arms carry simultaneous 2π sech pulses — the simulton regime where
coupled area-preserving solitons form.  Adding Doppler broadening via
thermal velocity classes forces the solver to run the full master-equation
at each velocity class and sum the result.

Features exercised in one solve:
  ✓ Large Hilbert space: 11×11 density matrix
  ✓ Two-field propagation (distinct field → coherence → field feedback loops)
  ✓ Hyperfine structure: Clebsch-Gordan-weighted coupling on both arms
  ✓ Two decay channels with per-channel factors (full Lindblad structure)
  ✓ Doppler velocity classes (most expensive scaling axis)
  ✓ Nonlinear (simulton) propagation → requires many ODE steps per z-slice

Run with:
    uv run pytest tests/bench_max.py --benchmark-only -v
"""

from maxwellbloch import hyperfine, mb_solve

# ---------------------------------------------------------------------------
# Build the hyperfine structure for the two arms of the V system.
#
# D1 arm: 5S₁/₂ F=1 (states 0-2) ↔ 5P₁/₂ F'=1 (states 3-5)
#   The Atom1e helper already assigns indices 0-2 to the ground F-level
#   and 3-5 to the excited F-level, so no remapping is needed.
#
# D2 arm: 5S₁/₂ F=1 (states 0-2) ↔ 5P₃/₂ F'=2 (states 6-10)
#   A separate Atom1e gives the ground states indices 0-2 and the excited
#   states indices 3-7.  We remap the excited indices by +3 to place them
#   at 6-10 in the combined 11-state Hilbert space.
# ---------------------------------------------------------------------------

# D1: 5S₁/₂ F=1 → 5P₁/₂ F'=1, σ⁺ (q=+1)
_d1 = hyperfine.Atom1e(element="Rb", isotope="87")
_d1.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))  # 5S₁/₂ F=1  → states 0-2
_d1.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))  # 5P₁/₂ F'=1 → states 3-5

D1_FIELD_CHANNELS = _d1.get_coupled_levels(F_level_idxs_a=(0,), F_level_idxs_b=(1,))
D1_FIELD_FACTORS = _d1.get_clebsch_hf_factors(
    F_level_idxs_a=(0,), F_level_idxs_b=(1,), q=1
).tolist()
D1_DECAY_CHANNELS = D1_FIELD_CHANNELS
D1_DECAY_FACTORS = _d1.get_decay_factors(
    F_level_idxs_a=(0,), F_level_idxs_b=(1,)
).tolist()

# D2: 5S₁/₂ F=1 → 5P₃/₂ F'=2, σ⁺ (q=+1)
_d2 = hyperfine.Atom1e(element="Rb", isotope="87")
_d2.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))  # 5S₁/₂ F=1  → states 0-2
_d2.add_F_level(hyperfine.LevelF(I=1.5, J=1.5, F=2))  # 5P₃/₂ F'=2 → states 3-7 (local)

_d2_coupled_raw = _d2.get_coupled_levels(F_level_idxs_a=(0,), F_level_idxs_b=(1,))
# Remap excited-state indices from local (3-7) → combined system (6-10)
D2_FIELD_CHANNELS = [[lo, hi + 3] for lo, hi in _d2_coupled_raw]
D2_FIELD_FACTORS = _d2.get_clebsch_hf_factors(
    F_level_idxs_a=(0,), F_level_idxs_b=(1,), q=1
).tolist()
D2_DECAY_CHANNELS = D2_FIELD_CHANNELS
D2_DECAY_FACTORS = _d2.get_decay_factors(
    F_level_idxs_a=(0,), F_level_idxs_b=(1,)
).tolist()

# Sech simulton pulse: area = rabi_freq × ampl × π × width = 2π → ampl = 2/width
_SECH_WIDTH = 0.3796628587572578
_SECH_AMPL = 2.0 / _SECH_WIDTH

# ---------------------------------------------------------------------------
# Full config
# ---------------------------------------------------------------------------

VEE_SIMULTON_RB87 = {
    "atom": {
        "num_states": 11,
        "decays": [
            {
                # D1: 5P₁/₂ F'=1 (states 3-5) → 5S₁/₂ F=1 (states 0-2)
                # Γ_D1 ≈ 5.75 MHz, normalised to 1.0
                "channels": D1_DECAY_CHANNELS,
                "rate": 1.0,
                "factors": D1_DECAY_FACTORS,
            },
            {
                # D2: 5P₃/₂ F'=2 (states 6-10) → 5S₁/₂ F=1 (states 0-2)
                # Γ_D2 ≈ 6.07 MHz ≈ 1.056 × Γ_D1; kept at 1.0 for simplicity
                "channels": D2_DECAY_CHANNELS,
                "rate": 1.0,
                "factors": D2_DECAY_FACTORS,
            },
        ],
        "fields": [
            {
                # D1 probe: 5S₁/₂ F=1 ↔ 5P₁/₂ F'=1, σ⁺
                "label": "d1_probe",
                "coupled_levels": D1_FIELD_CHANNELS,
                "factors": D1_FIELD_FACTORS,
                "detuning": 0.0,
                "detuning_positive": True,
                "rabi_freq": 1.0,
                "rabi_freq_t_func": "sech",
                "rabi_freq_t_args": {
                    "ampl": _SECH_AMPL,
                    "centre": 0.0,
                    "width": _SECH_WIDTH,
                },
            },
            {
                # D2 probe: 5S₁/₂ F=1 ↔ 5P₃/₂ F'=2, σ⁺
                "label": "d2_probe",
                "coupled_levels": D2_FIELD_CHANNELS,
                "factors": D2_FIELD_FACTORS,
                "detuning": 0.0,
                "detuning_positive": True,
                "rabi_freq": 1.0,
                "rabi_freq_t_func": "sech",
                "rabi_freq_t_args": {
                    "ampl": _SECH_AMPL,
                    "centre": 0.0,
                    "width": _SECH_WIDTH,
                },
            },
        ],
        # Equal population across 5S₁/₂ F=1 mF sublevels at t_min
        "initial_state": [1.0 / 3.0] * 3 + [0.0] * 8,
    },
    "t_min": -2.0,
    "t_max": 10.0,
    "t_steps": 100,
    "z_min": -0.2,
    "z_max": 1.2,
    "z_steps": 50,
    "z_steps_inner": 1,
    "interaction_strengths": [100.0, 100.0],
    "velocity_classes": {
        "thermal_width": 1.0,
        "thermal_delta_min": -3.0,
        "thermal_delta_max": 3.0,
        "thermal_delta_steps": 4,
        "thermal_delta_inner_min": 0.0,
        "thermal_delta_inner_max": 0.0,
        "thermal_delta_inner_steps": 0,
    },
}


def test_bench_max_vee_simulton_rb87(benchmark):
    """Rb 87 V-simulton: D1+D2 hyperfine, two 2π sech pulses, thermal.

    Hits every major code path in one solve:
    - 11×11 density matrix (both D1 and D2 hyperfine manifolds)
    - Two-field propagation with CG-weighted coupling on each arm
    - Two Lindblad decay channels with per-channel factors
    - Doppler velocity classes (runs mesolve independently per class)
    - Nonlinear (simulton) propagation regime
    """
    mbs = mb_solve.MBSolve(**VEE_SIMULTON_RB87)
    benchmark(mbs.mbsolve, recalc=True)
