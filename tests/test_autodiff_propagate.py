"""clerq.autodiff propagation vs ``PropagationSolver`` (the live oracle), and gradients."""

import copy
import glob
import json
import os

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import diffrax as dfx  # noqa: E402
import jax.numpy as jnp  # noqa: E402

from clerq import PropagationSolver  # noqa: E402
from clerq.autodiff._propagate import (  # noqa: E402
    propagate_fields,
    setup_from_solver,
)
from clerq.autodiff._pulses import pulse_function  # noqa: E402

pytestmark = pytest.mark.jax

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")
TIGHT = {"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000}
TOL = 1e-8
DEFAULT_INNER = 2  # PropagationSolver's default z_steps_inner


def _fixtures():
    out = {}
    for p in sorted(glob.glob(os.path.join(JSON_DIR, "*.json"))):
        with open(p) as f:
            d = json.load(f)
        if "z_steps" not in d:
            continue
        d = copy.deepcopy(d)
        # Cap only the time grid, to keep the oracle quick: reducing z_steps would
        # make the explicit z-scheme less stable than the fixture was designed to be.
        d["t_steps"] = min(d.get("t_steps", 100), 150)
        d.pop("savefile", None)
        out[os.path.basename(p)] = d
    return out


FIXTURES = _fixtures()

CUSTOM = {
    "system": {
        "num_states": 4,
        "energies": [0.0, 0.3, 1.2, 2.0],
        "decays": [
            {"rate": 0.7, "channels": [[0, 2], [1, 2]], "factors": [0.6, 0.8]},
            {"rate": 1.3, "channels": [[2, 3]]},
        ],
        "fields": [
            {
                "coupled_levels": [[0, 2], [1, 2]],
                "factors": [0.5, -0.8],
                "detuning": 0.4,
                "detuning_positive": False,
                "rabi_freq": 1.7,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.2, "centre": 1.0, "fwhm": 0.8},
            },
            {
                "coupled_levels": [[2, 3]],
                "detuning": -0.2,
                "factor_doppler_shift": 0.5,
                "rabi_freq": 0.9,
                "rabi_freq_t_func": "sech",
                "rabi_freq_t_args": {"ampl": 0.9, "centre": 1.5, "width": 0.4},
            },
        ],
        "initial_state": [0.5, 0.5, 0.0, 0.0],
    },
    "t_min": 0.0,
    "t_max": 3.0,
    "t_steps": 60,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_steps": 6,
    "z_steps_inner": 1,
    "interaction_strengths": [2.0, 1.0],
    "velocity_classes": {
        "thermal_delta_min": -1.0,
        "thermal_delta_max": 1.0,
        "thermal_delta_steps": 4,
        "thermal_width": 0.7,
    },
}


def _solver(d):
    ps = PropagationSolver.from_json_str(json.dumps(d))
    ps.build_opts(TIGHT)
    return ps


TIGHTER = {"atol": 1e-14, "rtol": 1e-12, "method": "adams", "nsteps": 1_000_000}


def _clerq_omegas(d, step, opts=TIGHT):
    ps = PropagationSolver.from_json_str(json.dumps(d))
    ps.build_opts(opts)
    ps.solve(step=step, recalc=True, progress=False)
    return ps._Omegas_zt


def _rel(a, b):
    return float(np.abs(a - b).max() / np.abs(b).max())


def _reference(d, step, opts=TIGHT):
    """clerq's own result for the same problem and scheme.

    Adams-Bashforth with ``z_steps_inner > 1`` stored every slice one inner step early
    (and read the density one inner step ahead) until clerq #328; with that fix,
    ``PropagationSolver`` and ``to_autodiff`` agree for every ``z_steps_inner``.
    """
    return _clerq_omegas(d, step, opts)


_CONDITIONED = {}
_NATIVE = {}  # (key, step) -> (problem, clerq at TIGHT, clerq at TIGHTER) at scale 1


def _conditioned(key, d, step):
    """A problem whose clerq solution does not depend on its own tolerance.

    Some fixtures are used in clerq's tests at interaction strengths where the
    explicit z-scheme is exponentially unstable (clerq differs from *itself* by
    percent between rtol 1e-10 and 1e-12), so no two solvers can agree there. Return
    ``(problem, reference, scale)`` with the interaction strengths reduced by
    ``scale`` until clerq agrees with itself to 1e-9 at two tolerances.
    """
    if (key, step) not in _CONDITIONED:
        result = None
        for scale in (1.0, 0.3, 0.1, 0.03, 0.01, 0.003):
            ds = {
                **d,
                "interaction_strengths": [
                    s * scale for s in d["interaction_strengths"]
                ],
            }
            a = _reference(ds, step, TIGHT)
            b = _reference(ds, step, TIGHTER)
            if scale == 1.0:
                _NATIVE[(key, step)] = (ds, a, b)
            if _rel(a, b) < 1e-9:
                result = (ds, a, scale)
                break
        _CONDITIONED[(key, step)] = result
    result = _CONDITIONED[(key, step)]
    if result is None:
        pytest.skip("no interaction strength gives a well-conditioned problem")
    return result


def _has_medium(d):
    """False for fixtures with zero number density (the field cannot change)."""
    return d.get("num_density_z_args", {"ampl": 1.0}).get("ampl", 1.0) != 0.0


def test_oracle_set_is_not_vacuous():
    live = {k: d for k, d in FIXTURES.items() if _has_medium(d)}
    assert len(FIXTURES) >= 10
    assert len(live) >= 8, "too few fixtures with an actual medium"
    inners = {d.get("z_steps_inner", DEFAULT_INNER) for d in live.values()}
    assert 1 in inners and max(inners) > 1
    assert any(len(d["system"]["fields"]) >= 2 for d in live.values())
    assert any("velocity_classes" in d for d in live.values())


@pytest.mark.parametrize("step", ["euler", "ab"])
@pytest.mark.parametrize("name", list(FIXTURES))
def test_matches_propagation_solver(name, step):
    try:
        d, ref, scale = _conditioned(name, FIXTURES[name], step)
    except pytest.skip.Exception:
        raise
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"clerq cannot solve this fixture itself: {type(e).__name__}")
    setup = setup_from_solver(_solver(d))
    got = np.asarray(propagate_fields(setup, step=step, rtol=1e-11, atol=1e-13))
    assert got.shape == ref.shape
    assert _rel(got, ref) < TOL, (scale, _rel(got, ref))
    # non-vacuity: where there is a medium it must actually change the field
    if _has_medium(d):
        assert np.abs(ref[:, -1] - ref[:, 0]).max() > 1e-8 * np.abs(ref[:, 0]).max()
    else:
        assert np.array_equal(got[:, -1], got[:, 0])  # no medium: field unchanged


@pytest.mark.parametrize("step", ["euler", "ab"])
@pytest.mark.parametrize("name", list(FIXTURES))
def test_matches_at_native_coupling_within_the_oracles_noise(name, step):
    """At the fixture's own interaction strength the explicit z-scheme can amplify
    solver tolerance, so clerq differs from itself between two tolerances. Our result
    must agree with clerq's tighter run to within that self-disagreement; a wrong
    convention is an O(1) error, far outside it."""
    try:
        _conditioned(name, FIXTURES[name], step)  # populates _NATIVE
    except pytest.skip.Exception:
        pass
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"clerq cannot solve this fixture itself: {type(e).__name__}")
    d, a, b = _NATIVE[(name, step)]
    noise = _rel(a, b)
    setup = setup_from_solver(_solver(d))
    got = np.asarray(propagate_fields(setup, step=step, rtol=1e-12, atol=1e-14))
    err = _rel(got, b)
    assert err <= max(TOL, 3.0 * noise), (err, noise)
    if _has_medium(d):  # and it must be a meaningful comparison
        assert np.abs(b[:, -1] - b[:, 0]).max() > 1e-8 * np.abs(b[:, 0]).max()
        assert err < 0.1, "the agreement bound is too loose to mean anything"


@pytest.mark.parametrize(
    "step,inner", [("euler", 1), ("euler", 2), ("ab", 1), ("ab", 2)]
)
def test_multi_field_velocity_classes_energies_doppler_factor(step, inner):
    d, ref, scale = _conditioned(
        f"custom-inner{inner}", {**CUSTOM, "z_steps_inner": inner}, step
    )
    assert scale == 1.0  # designed to be stable: no reduction of the coupling
    setup = setup_from_solver(_solver(d))
    got = np.asarray(propagate_fields(setup, step=step, rtol=1e-11, atol=1e-13))
    assert _rel(got, ref) < TOL, _rel(got, ref)
    assert np.abs(ref[:, -1] - ref[:, 0]).max() > 1e-3  # a real effect
    assert setup.ops.num_fields == 2 and len(setup.shifts) == 5


def test_slices_are_at_their_nominal_z_with_inner_steps():
    """With inner steps, the stored slice j is the same scheme as slice j*inner of a
    run that stores every inner step (the property clerq #328 restored): a problem
    with z_steps=N, inner=2 equals the z_steps=2N, inner=1 run sampled every 2nd."""
    d2 = {**CUSTOM, "z_steps": 3, "z_steps_inner": 2}
    d1 = {**CUSTOM, "z_steps": 6, "z_steps_inner": 1}
    setup2 = setup_from_solver(_solver(d2))
    setup1 = setup_from_solver(_solver(d1))
    a = np.asarray(propagate_fields(setup2, step="ab", rtol=1e-11, atol=1e-13))
    b = np.asarray(propagate_fields(setup1, step="ab", rtol=1e-11, atol=1e-13))
    assert _rel(a, b[:, ::2, :]) < TOL, _rel(a, b[:, ::2, :])
    # and clerq itself agrees with the same property
    ca = _clerq_omegas(d2, "ab")
    cb = _clerq_omegas(d1, "ab")
    assert _rel(ca, cb[:, ::2, :]) < TOL, _rel(ca, cb[:, ::2, :])


# ---- guards ---------------------------------------------------------------------


def test_unsupported_pulse_is_rejected():
    with pytest.raises(NotImplementedError, match="sinc"):
        pulse_function("sinc", {})


def test_unknown_step_and_strategy_are_rejected():
    setup = setup_from_solver(_solver(CUSTOM))
    with pytest.raises(ValueError, match="step"):
        propagate_fields(setup, step="rk4")
    with pytest.raises(ValueError, match="strategy"):
        propagate_fields(setup, strategy="s9")


# ---- gradients --------------------------------------------------------------------

SMALL = {**CUSTOM, "t_steps": 40, "z_steps": 4}
N_FIXED = 2003  # not a multiple of the 40 output intervals


def _problem():
    setup = setup_from_solver(_solver(SMALL))
    n_f, n_t = setup.ops.num_fields, len(setup.tlist)
    omega_in = np.asarray(propagate_fields(setup, rtol=1e-10, atol=1e-12)[:, 0, :])
    rng = np.random.default_rng(5)
    probe = jnp.asarray(rng.normal(size=(n_f, n_t)))
    wt = jnp.asarray(np.exp(-0.5 * ((setup.tlist - 1.5) / 0.8) ** 2))

    def unpack(x):
        o = x[: n_f * n_t].reshape(n_f, n_t)
        g = x[n_f * n_t : n_f * n_t + n_f]
        det = x[n_f * n_t + n_f : n_f * n_t + 2 * n_f]
        scale = x[-1]
        return o, g, det, scale

    def loss(x, **kw):
        o, g, det, scale = unpack(x)
        out = propagate_fields(
            setup, omega0=o, g=g, detunings=det, decay_scale=scale, **kw
        )
        z_last = out[:, -1, :]
        return jnp.sum(jnp.abs(z_last) ** 2 * wt) + jnp.sum(jnp.imag(z_last) * probe)

    x0 = np.concatenate(
        [omega_in.real.ravel(), 2.0 * np.array([1.0, 0.5]), setup.ops.detunings, [1.2]]
    )
    idx = []
    for k in range(n_f):
        base = k * n_t
        idx += [base + i for i in (0, 1, 10, 20, 30, n_t - 2, n_t - 1)]
    idx += [n_f * n_t + k for k in range(n_f)]  # g
    idx += [n_f * n_t + n_f + k for k in range(n_f)]  # detunings
    idx += [n_f * n_t + 2 * n_f]  # decay scale
    return loss, jnp.asarray(x0), idx


def _central_fd(f, x, idx, h):
    x = np.asarray(x, dtype=float)
    out = np.empty(len(idx))
    for j, i in enumerate(idx):
        e = np.zeros_like(x)
        e[i] = h
        out[j] = (float(f(jnp.asarray(x + e))) - float(f(jnp.asarray(x - e)))) / (2 * h)
    return out


def _grel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def test_gradients_match_finite_differences():
    loss, x0, idx = _problem()
    fixed = lambda x: loss(x, n_steps=N_FIXED)  # noqa: E731
    g_fixed = np.asarray(jax.jit(jax.grad(fixed))(x0))
    fd = _central_fd(jax.jit(fixed), x0, idx, 1e-4)
    assert _grel(g_fixed[idx], fd) < 1e-5, _grel(g_fixed[idx], fd)
    assert np.abs(g_fixed[idx]).max() > 1e-4  # non-vacuous

    adaptive = lambda x: loss(x, rtol=1e-11, atol=1e-13)  # noqa: E731
    g_pid = np.asarray(jax.jit(jax.grad(adaptive))(x0))
    assert _grel(g_pid, g_fixed) < 1e-4, _grel(g_pid, g_fixed)


def test_strategies_give_the_same_gradient():
    loss, x0, _ = _problem()
    grads = {}
    for s, kw in [("s2", {}), ("s3", {}), ("s4", {"checkpoints": 16})]:
        f = lambda x, s=s, kw=kw: loss(  # noqa: E731
            x, strategy=s, rtol=1e-10, atol=1e-12, **kw
        )
        grads[s] = np.asarray(jax.jit(jax.grad(f))(x0))
    assert _grel(grads["s3"], grads["s2"]) < 1e-9
    assert _grel(grads["s4"], grads["s2"]) < 1e-9


def test_reverse_mode_matches_forward_mode():
    loss, x0, _ = _problem()
    rev = lambda x: loss(x, rtol=1e-10, atol=1e-12)  # noqa: E731
    fwd = lambda x: loss(  # noqa: E731
        x, rtol=1e-10, atol=1e-12, adjoint=dfx.ForwardMode()
    )
    g_rev = np.asarray(jax.jit(jax.grad(rev))(x0))
    g_fwd = np.asarray(jax.jit(jax.jacfwd(fwd))(x0))
    assert _grel(g_rev, g_fwd) < 1e-8, _grel(g_rev, g_fwd)
