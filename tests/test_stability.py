import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from clerq import PropagationSolver
from clerq.exceptions import PropagationInstabilityError

JSON = Path(__file__).parent / "json" / "mbs_two_sech_2pi.json"


def make(strength=10.0, z_steps_inner=2, z_steps=10):
    d = json.loads(JSON.read_text())
    d.update(
        savefile=None,
        z_steps=z_steps,
        z_steps_inner=z_steps_inner,
        interaction_strengths=[strength],
    )
    return PropagationSolver(**d)


@pytest.mark.parametrize("step", ["euler", "ab"])
def test_coarse_step_warns(step):
    s = make(z_steps_inner=2)
    with pytest.warns(UserWarning, match="z-step may be too coarse"):
        s.solve(step=step, progress=False)
    assert s.step_error_estimate > 0.1


def test_fine_step_is_quiet():
    s = make(z_steps_inner=32)
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
    s = make(z_steps_inner="auto")
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
