import json
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from clerq import PropagationSolver
from clerq.exceptions import PropagationInstabilityError

JSON = Path(__file__).parent / "json" / "mbs_two_sech_2pi.json"


def make(strength=10.0, z_steps_inner=2, z_steps=10, **overrides):
    d = json.loads(JSON.read_text())
    d.update(
        savefile=None,
        z_steps=z_steps,
        z_steps_inner=z_steps_inner,
        interaction_strengths=[strength],
    )
    d.update(overrides)
    return PropagationSolver(**d)


# The pulse is delayed by the medium; a wide window keeps it from being clipped.
WIDE = dict(t_max=20.0, t_steps=240)


@pytest.mark.parametrize("step", ["euler", "ab"])
def test_coarse_step_warns(step):
    s = make(z_steps_inner=2)
    with pytest.warns(UserWarning, match="z-step may be too coarse"):
        s.solve(step=step, progress=False)
    assert s.step_error_estimate > 0.1


def test_fine_step_is_quiet():
    s = make(z_steps_inner=32, **WIDE)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Omegas, _ = s.solve(step="ab", progress=False)
    # A 2pi sech pulse propagates undistorted: peak stays near its input.
    assert np.abs(Omegas[0, -1]).max() == pytest.approx(
        np.abs(Omegas[0, 0]).max(), rel=0.02
    )


def test_blowup_stops_with_error():
    s = make()
    with pytest.raises(PropagationInstabilityError) as exc:
        s.solve(progress=False, blowup_factor=2.0)
    assert exc.value.suggested_z_steps_inner > s.z_steps_inner
    assert exc.value.z is not None


def test_integrator_failure_is_wrapped():
    s = make(strength=3000.0)
    with pytest.raises(PropagationInstabilityError, match="integrator failed"):
        s.solve(progress=False)


def test_check_stability_false_suppresses():
    s = make()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        s.solve(progress=False, check_stability=False, blowup_factor=2.0)


def test_auto_z_steps_inner():
    s = make(z_steps_inner="auto", **WIDE)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Omegas, _ = s.solve(progress=False)
    assert s.z_steps_inner > 2
    assert s.step_error_estimate <= 0.1
    assert np.abs(Omegas[0, -1]).max() == pytest.approx(
        np.abs(Omegas[0, 0]).max(), rel=0.1
    )


def test_auto_gives_up_at_max():
    s = make(z_steps_inner="auto")
    with pytest.raises(PropagationInstabilityError, match="max_z_steps_inner"):
        s.solve(progress=False, max_z_steps_inner=4)


class TestTimeResolution:
    def test_decimation_error_converges_as_fourth_power(self):
        def err(n):
            t = np.linspace(-4, 4, n + 1)
            y = np.exp(-(t**2)) * np.cos(3 * t)
            return PropagationSolver._decimation_error(y, np.array(1.0))

        ratio = err(80) / err(160)
        assert 12 < ratio < 24

    def test_decimation_error_needs_enough_samples(self):
        assert PropagationSolver._decimation_error(np.ones((2, 5)), np.ones(2)) == 0.0

    def test_coarse_time_grid_warns(self):
        s = make(z_steps_inner=32, **{**WIDE, "t_steps": 60})
        with pytest.warns(UserWarning, match="time grid may be too coarse"):
            s.solve(progress=False)
        assert s.t_step_error_estimate > 0.1

    def test_error_estimate_falls_with_t_steps(self):
        errs = []
        for t_steps in (120, 240):
            s = make(z_steps_inner=32, **{**WIDE, "t_steps": t_steps})
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s.solve(progress=False)
            errs.append(s.t_step_error_estimate)
        assert 8 < errs[0] / errs[1] < 40

    def test_following_the_suggestion_converges(self):
        # Far from the asymptotic regime the suggestion is a lower bound, so
        # following it a few times is what clears the warning.
        t_steps = 60
        for _ in range(4):
            s = make(z_steps_inner=32, **{**WIDE, "t_steps": t_steps})
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                s.solve(progress=False)
            msgs = [str(w.message) for w in rec if "time grid" in str(w.message)]
            if not msgs:
                break
            nxt = int(re.search(r"t_steps=(\d+)", msgs[0]).group(1))
            assert nxt > t_steps
            t_steps = nxt
        assert not msgs
        assert t_steps > 60

    def test_clipped_window_warns(self):
        # t_max=10 cuts off the delayed pulse.
        s = make(z_steps_inner=32, t_steps=240)
        with pytest.warns(UserWarning, match="edge of the time window"):
            s.solve(progress=False)
        assert s.window_edge_estimate > 0.2

    def test_wide_window_is_quiet(self):
        s = make(z_steps_inner=32, **WIDE)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            s.solve(progress=False)
        assert s.window_edge_estimate < 0.2
        assert s.t_step_error_estimate < 0.1

    def test_check_stability_false_skips_time_checks(self):
        s = make(z_steps_inner=32, t_steps=40)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            s.solve(progress=False, check_stability=False)
        assert s.t_step_error_estimate == 0.0
        assert s.window_edge_estimate == 0.0

    def test_suggest_t_steps_resolves_the_input_pulse(self):
        s = make(**WIDE)
        n = s.suggest_t_steps()
        # Sech pulse of width 0.38 over a window of 22: tens to a few hundred.
        assert 50 < n < 400

    def test_suggest_t_steps_grows_with_detuning(self):
        base = make(**WIDE).suggest_t_steps()
        d = json.loads(JSON.read_text())
        d["system"]["fields"][0]["detuning"] = 20.0
        d.update(savefile=None, interaction_strengths=[10.0], **WIDE)
        detuned = PropagationSolver(**d).suggest_t_steps()
        # 10 points per period of 20 over a window of 22 is at least ~4400 steps.
        assert detuned >= 4000
        assert detuned > 10 * base
