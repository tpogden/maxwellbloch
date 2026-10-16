"""Clear errors for problems and options the autodiff backend cannot handle."""

import json

import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

from clerq import PropagationSolver  # noqa: E402

pytestmark = pytest.mark.jax


def _problem(fields=True):
    system = {
        "num_states": 2,
        "decays": [{"rate": 1.0, "channels": [[0, 1]]}],
        "fields": [],
    }
    if fields:
        system["fields"] = [
            {
                "coupled_levels": [[0, 1]],
                "rabi_freq": 1.0,
                "detuning": 0.0,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.0, "centre": 1.0, "fwhm": 1.0},
            }
        ]
    return {
        "system": system,
        "t_min": 0.0,
        "t_max": 3.0,
        "t_steps": 20,
        "z_steps": 2,
        "z_steps_inner": 1,
        "interaction_strengths": [1.0] if fields else [],
    }


def _solver(**kw):
    return PropagationSolver.from_json_str(json.dumps(_problem(**kw)))


def test_system_without_fields_is_rejected_clearly():
    with pytest.raises(ValueError, match="at least one coupled field"):
        _solver(fields=False).to_autodiff()


def test_s4_requires_checkpoints():
    ad = _solver().to_autodiff(strategy="s4")
    with pytest.raises(ValueError, match="needs `checkpoints`"):
        ad.propagate(ad.default_params())


@pytest.mark.parametrize("bad", [0, -3, 2.5])
def test_s4_rejects_a_non_positive_or_fractional_checkpoint_count(bad):
    ad = _solver().to_autodiff(strategy="s4", checkpoints=bad)
    with pytest.raises(ValueError, match="positive integer"):
        ad.propagate(ad.default_params())


@pytest.mark.parametrize("strategy", ["s2", "s3"])
def test_checkpoints_is_rejected_for_the_other_strategies(strategy):
    ad = _solver().to_autodiff(strategy=strategy, checkpoints=8)
    with pytest.raises(ValueError, match="only applies to strategy='s4'"):
        ad.propagate(ad.default_params())


def test_s4_with_checkpoints_runs():
    ad = _solver().to_autodiff(strategy="s4", checkpoints=8)
    out = ad.propagate(ad.default_params())
    assert out.shape[1] == 3
