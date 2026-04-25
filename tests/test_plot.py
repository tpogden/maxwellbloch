"""Smoke tests for maxwellbloch.plot primitives.

Tests verify that every public primitive:
  - returns a plotly Figure
  - has the expected trace type
  - has layout attributes set (title, axis labels)

No visual correctness is checked — that belongs in the usage notebook.
"""

import unittest

import pytest

pytest.importorskip("plotly", reason="plotly not installed; install maxwellbloch[plot]")

import numpy as np
import plotly.graph_objects as go

from maxwellbloch import mb_solve, plot

# Minimal two-level two-field config used by all field/state tests.
_TWO_LEVEL_JSON = """
{
  "atom": {
    "num_states": 2,
    "decays": [{"channels": [[0, 1]], "rate": 1.0}],
    "fields": [
      {
        "label": "probe",
        "coupled_levels": [[0, 1]],
        "rabi_freq": 0.01,
        "rabi_freq_t_func": "gaussian",
        "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": 1.0}
      }
    ]
  },
  "t_min": -2.0,
  "t_max": 4.0,
  "t_steps": 60,
  "z_min": 0.0,
  "z_max": 1.0,
  "z_steps": 5,
  "z_steps_inner": 2,
  "interaction_strengths": [0.1]
}
"""

# Three-level ladder — two fields, for spectrum_overlay test.
_THREE_LEVEL_JSON = """
{
  "atom": {
    "num_states": 3,
    "decays": [{"channels": [[0, 1], [1, 2]], "rate": 1.0}],
    "fields": [
      {
        "label": "probe",
        "coupled_levels": [[0, 1]],
        "rabi_freq": 0.01,
        "rabi_freq_t_func": "gaussian",
        "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": 1.0}
      },
      {
        "label": "coupling",
        "coupled_levels": [[1, 2]],
        "rabi_freq": 5.0,
        "rabi_freq_t_func": "ramp_onoff",
        "rabi_freq_t_args": {"ampl": 1.0, "fwhm": 0.2, "on": -2.0, "off": 4.0}
      }
    ]
  },
  "t_min": -2.0,
  "t_max": 4.0,
  "t_steps": 60,
  "z_min": 0.0,
  "z_max": 1.0,
  "z_steps": 5,
  "z_steps_inner": 2,
  "interaction_strengths": [0.1, 0.1]
}
"""


def _solve_two_level():
    mbs = mb_solve.MBSolve.from_json_str(_TWO_LEVEL_JSON)
    mbs.mbsolve(recalc=True)
    return mbs


def _solve_three_level():
    mbs = mb_solve.MBSolve.from_json_str(_THREE_LEVEL_JSON)
    mbs.mbsolve(recalc=True)
    return mbs


# Shared solved instances (module-level to avoid re-solving per test).
_MBS2 = None
_MBS3 = None


def _get_mbs2():
    global _MBS2
    if _MBS2 is None:
        _MBS2 = _solve_two_level()
    return _MBS2


def _get_mbs3():
    global _MBS3
    if _MBS3 is None:
        _MBS3 = _solve_three_level()
    return _MBS3


class TestFieldSpacetime(unittest.TestCase):
    def test_returns_figure(self):
        fig = plot.field_spacetime(_get_mbs2())
        self.assertIsInstance(fig, go.Figure)

    def test_has_heatmap_trace(self):
        fig = plot.field_spacetime(_get_mbs2())
        self.assertIsInstance(fig.data[0], go.Heatmap)

    def test_heatmap_shape(self):
        mbs = _get_mbs2()
        fig = plot.field_spacetime(mbs)
        z_vals = fig.data[0].z
        self.assertEqual(np.array(z_vals).shape[0], len(mbs.zlist))

    def test_title_contains_field_label(self):
        fig = plot.field_spacetime(_get_mbs2())
        self.assertIn("probe", fig.layout.title.text)

    def test_axes_labelled(self):
        fig = plot.field_spacetime(_get_mbs2())
        self.assertIn("t", fig.layout.xaxis.title.text)
        self.assertIn("z", fig.layout.yaxis.title.text)

    def test_custom_colorscale(self):
        fig_v = plot.field_spacetime(_get_mbs2(), colorscale="Viridis")
        fig_c = plot.field_spacetime(_get_mbs2(), colorscale="Cividis")
        # Plotly resolves names to tuples; they should differ between scales
        self.assertNotEqual(fig_v.data[0].colorscale, fig_c.data[0].colorscale)


class TestFieldEnvelope(unittest.TestCase):
    def test_returns_figure(self):
        fig = plot.field_envelope(_get_mbs2())
        self.assertIsInstance(fig, go.Figure)

    def test_default_two_traces(self):
        fig = plot.field_envelope(_get_mbs2())
        self.assertEqual(len(fig.data), 2)

    def test_single_z_index(self):
        fig = plot.field_envelope(_get_mbs2(), z_indices=0)
        self.assertEqual(len(fig.data), 1)

    def test_trace_length_matches_tlist(self):
        mbs = _get_mbs2()
        fig = plot.field_envelope(mbs, z_indices=0)
        self.assertEqual(len(fig.data[0].x), len(mbs.tlist))


class TestPulseArea(unittest.TestCase):
    def test_returns_figure(self):
        fig = plot.pulse_area(_get_mbs2())
        self.assertIsInstance(fig, go.Figure)

    def test_has_scatter_trace(self):
        fig = plot.pulse_area(_get_mbs2())
        self.assertIsInstance(fig.data[0], go.Scatter)

    def test_trace_length_matches_zlist(self):
        mbs = _get_mbs2()
        fig = plot.pulse_area(mbs)
        self.assertEqual(len(fig.data[0].x), len(mbs.zlist))

    def test_area_decreases_for_absorber(self):
        mbs = _get_mbs2()
        fig = plot.pulse_area(mbs)
        areas = fig.data[0].y
        # Probe loses energy to absorber — area should decrease z_min → z_max
        self.assertGreater(areas[0], areas[-1])


class TestSpectrum(unittest.TestCase):
    def test_returns_figure(self):
        fig = plot.spectrum(_get_mbs2())
        self.assertIsInstance(fig, go.Figure)

    def test_has_scatter_trace(self):
        fig = plot.spectrum(_get_mbs2())
        self.assertIsInstance(fig.data[0], go.Scatter)

    def test_with_dispersion(self):
        fig = plot.spectrum(_get_mbs2(), show_dispersion=True)
        self.assertEqual(len(fig.data), 2)

    def test_freq_axis_length(self):
        mbs = _get_mbs2()
        fig = plot.spectrum(mbs)
        self.assertEqual(len(fig.data[0].x), len(mbs.tlist))


class TestSpectrumOverlay(unittest.TestCase):
    def test_returns_figure(self):
        mbs = _get_mbs2()
        fig = plot.spectrum_overlay([mbs, mbs], labels=["a", "b"])
        self.assertIsInstance(fig, go.Figure)

    def test_trace_count(self):
        mbs = _get_mbs2()
        fig = plot.spectrum_overlay([mbs, mbs])
        self.assertEqual(len(fig.data), 2)


class TestPopulation(unittest.TestCase):
    def test_returns_figure(self):
        fig = plot.population(_get_mbs2())
        self.assertIsInstance(fig, go.Figure)

    def test_default_all_states(self):
        mbs = _get_mbs2()
        fig = plot.population(mbs)
        self.assertEqual(len(fig.data), mbs.atom.num_states)

    def test_single_state(self):
        fig = plot.population(_get_mbs2(), state_indices=0)
        self.assertEqual(len(fig.data), 1)

    def test_populations_sum_to_one(self):
        mbs = _get_mbs2()
        fig = plot.population(mbs)
        total = sum(np.array(t.y) for t in fig.data)
        np.testing.assert_allclose(total, 1.0, atol=1e-6)


class TestPopulationSpacetime(unittest.TestCase):
    def test_returns_figure(self):
        fig = plot.population_spacetime(_get_mbs2(), state_idx=0)
        self.assertIsInstance(fig, go.Figure)

    def test_has_heatmap_trace(self):
        fig = plot.population_spacetime(_get_mbs2(), state_idx=0)
        self.assertIsInstance(fig.data[0], go.Heatmap)

    def test_colorbar_range_default_auto(self):
        fig = plot.population_spacetime(_get_mbs2(), state_idx=0)
        self.assertIsNone(fig.data[0].zmin)
        self.assertIsNone(fig.data[0].zmax)

    def test_colorbar_range_explicit(self):
        fig = plot.population_spacetime(_get_mbs2(), state_idx=0, zmin=0.0, zmax=1.0)
        self.assertEqual(fig.data[0].zmin, 0.0)
        self.assertEqual(fig.data[0].zmax, 1.0)


class TestCoherence(unittest.TestCase):
    def test_returns_figure(self):
        fig = plot.coherence(_get_mbs2(), 0, 1)
        self.assertIsInstance(fig, go.Figure)

    def test_abs_component(self):
        mbs = _get_mbs2()
        fig = plot.coherence(mbs, 0, 1, component="abs")
        self.assertTrue(np.all(np.array(fig.data[0].y) >= 0))

    def test_invalid_component_raises(self):
        with self.assertRaises(ValueError):
            plot.coherence(_get_mbs2(), 0, 1, component="phase")


if __name__ == "__main__":
    unittest.main()
