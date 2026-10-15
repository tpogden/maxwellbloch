"""clerq.autodiff single-slice solver vs ``LindbladSolver`` / clerq's own drive, and gradients.

clerq is the live oracle: the same problem definition is solved by clerq and by the
autodiff solver and the density matrices are compared.
"""

import copy
import dataclasses
import glob
import json
import math
import os
import zlib

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import diffrax as dfx  # noqa: E402
import jax.numpy as jnp  # noqa: E402

from clerq import LindbladSolver, PropagationSolver  # noqa: E402
from clerq.autodiff._drive import CubicDrive  # noqa: E402
from clerq.autodiff._lindblad import solve_slice  # noqa: E402
from clerq.autodiff._operators import build_operators  # noqa: E402
from clerq.system import System  # noqa: E402

pytestmark = pytest.mark.jax

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")
# clerq's own solver error is what limits the comparison: at rtol 1e-10 it is ~5e-9
# (and ~1e-8 on some platforms), at rtol 1e-12 it is ~2e-10, so the 1e-8 bar below
# leaves a comfortable margin instead of a marginal one.
TIGHT = {"atol": 1e-14, "rtol": 1e-12, "method": "adams", "nsteps": 1_000_000}
TOL = 1e-8


def _load(path):
    with open(path) as f:
        d = json.load(f)
    sd = d["system"] if "system" in d else (d if "num_states" in d else None)
    if sd is None or not sd.get("fields"):
        return None
    grid = (d.get("t_min", 0.0), d.get("t_max", 1.0), max(d.get("t_steps", 100), 100))
    return sd, grid


FIXTURES = {
    os.path.basename(p): v
    for p in sorted(glob.glob(os.path.join(JSON_DIR, "*.json")))
    if (v := _load(p)) is not None
}

CUSTOM = (
    {
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
                "rabi_freq": 0.9,
                "rabi_freq_t_func": "sech",
                "rabi_freq_t_args": {"ampl": 0.9, "centre": 1.5, "width": 0.4},
            },
        ],
        "initial_state": [0.5, 0.5, 0.0, 0.0],
    },
    (0.0, 3.0, 90),
)
ALL = {**FIXTURES, "custom": CUSTOM}


def _jax_pulse(field):
    """JAX version of clerq's gaussian / sech pulse for a field dict (cyclic units)."""
    name = field.get("rabi_freq_t_func")
    args = field.get("rabi_freq_t_args", {})
    r0 = field.get("rabi_freq", 1.0)
    if name == "gaussian":
        fwhm, centre = args["fwhm"], args["centre"]
        if "ampl" in args:
            ampl = args["ampl"]
        else:
            ampl = (
                args["n_pi"]
                * math.sqrt(4 * math.pi * math.log(2) / fwhm**2)
                / (2 * math.pi)
            )
        k = 4 * math.log(2)
        return lambda t: r0 * ampl * jnp.exp(-k * ((t - centre) / fwhm) ** 2)
    if name == "sech":
        centre = args["centre"]
        width = args["width"] if "width" in args else args["fwhm"] / 2.6339157938
        ampl = args["ampl"] if "ampl" in args else args["n_pi"] / width / (2 * math.pi)
        return lambda t: (
            r0
            * ampl
            * 2.0
            / (jnp.exp((t - centre) / width) + jnp.exp(-(t - centre) / width))
        )
    return None


def _is_analytic(sd):
    return all(_jax_pulse(f) is not None for f in sd["fields"])


ANALYTIC = {k: v for k, v in ALL.items() if _is_analytic(v[0])}


def _samples(n_fields, tlist, seed):
    """Smooth, non-Gaussian cyclic Rabi-frequency samples, one row per field."""
    rng = np.random.default_rng(seed)
    t0, t1 = tlist[0], tlist[-1]
    out = []
    for k in range(n_fields):
        u = (tlist - t0) / (t1 - t0)
        centre = 0.3 + 0.4 * (k % 2)
        base = 0.9 * np.exp(-(((u - centre) / 0.18) ** 2)) + 0.25 * np.sin(5.0 * u + k)
        out.append(
            base * np.exp(-0.5 * ((u - 0.5) / 0.6) ** 2)
            + 0.02 * rng.normal(size=u.shape)
        )
    return np.array(out)


def _clerq_lindblad(sd, grid):
    s = LindbladSolver(
        system=copy.deepcopy(sd),
        t_min=grid[0],
        t_max=grid[1],
        t_steps=grid[2],
        opts=TIGHT,
    )
    s.solve(save=False)
    return s.tlist, s.states_t()


def _clerq_interpolated(sd, grid, samples_cyclic):
    """Solve with the field given as samples, interpolated exactly as clerq does."""
    ps = PropagationSolver(
        system=copy.deepcopy(sd),
        t_min=grid[0],
        t_max=grid[1],
        t_steps=grid[2],
        z_steps=1,
        opts=TIGHT,
    )
    omegas = (2 * np.pi * samples_cyclic).astype(complex)  # clerq stores angular
    ps.system.H_Omega_list = ps._build_intp_H_Omega_list(omegas)
    ps.system.mesolve(ps.tlist, options=ps.opts, recalc=True)
    return ps.tlist, ps.system.states_t()


def test_oracle_sets_are_not_vacuous():
    assert len(FIXTURES) >= 15
    assert len(ANALYTIC) >= 8
    assert any(len(v[0]["fields"]) >= 2 for v in ANALYTIC.values()), (
        "no analytic multi-field"
    )
    assert any(len(v[0]["fields"]) >= 2 for v in FIXTURES.values()), (
        "no sampled multi-field"
    )
    assert any(v[0].get("energies") for v in ALL.values()), "no bare energies"


@pytest.mark.parametrize("name", list(ANALYTIC))
def test_analytic_drive_matches_lindblad_solver(name):
    sd, grid = ANALYTIC[name]
    tlist, ref = _clerq_lindblad(sd, grid)
    ops = build_operators(System(**copy.deepcopy(sd)))
    pulses = [_jax_pulse(f) for f in sd["fields"]]
    rabi = lambda t: jnp.stack([p(t) for p in pulses])  # noqa: E731
    rho = solve_slice(ops, rabi, ops.detunings, tlist, rtol=1e-11, atol=1e-13)
    got = np.asarray(rho)
    assert np.abs(got - ref).max() < TOL, np.abs(got - ref).max()
    # non-vacuity: the drive must actually move the state
    assert np.abs(got - got[0]).max() > 1e-4


@pytest.mark.parametrize("name", list(ALL))
def test_sampled_drive_matches_clerq_interpolation(name):
    sd, grid = ALL[name]
    n_f = len(sd["fields"])
    tlist = np.linspace(grid[0], grid[1], grid[2] + 1)
    # deterministic: str.__hash__ is randomised per process, which made this flaky
    samples = _samples(n_f, tlist, seed=zlib.crc32(name.encode()) % 1000)
    t_ref, ref = _clerq_interpolated(sd, grid, samples)
    ops = build_operators(System(**copy.deepcopy(sd)))
    rabi = CubicDrive(t_ref)(jnp.asarray(samples))
    got = np.asarray(
        solve_slice(ops, rabi, ops.detunings, t_ref, rtol=1e-11, atol=1e-13)
    )
    assert np.abs(got - ref).max() < TOL, np.abs(got - ref).max()
    assert np.abs(got - got[0]).max() > 1e-4


def test_sparse_and_dense_dissipators_agree():
    sd, grid = CUSTOM
    ops = build_operators(System(**copy.deepcopy(sd)))
    assert ops.sparse is not None
    dense = dataclasses.replace(ops, sparse=None)
    tlist = np.linspace(grid[0], grid[1], 61)
    rabi = lambda t: jnp.array([0.7, 0.4]) * jnp.exp(-((t - 1.2) ** 2))  # noqa: E731
    a = solve_slice(ops, rabi, ops.detunings, tlist, rtol=1e-11, atol=1e-13)
    b = solve_slice(dense, rabi, ops.detunings, tlist, rtol=1e-11, atol=1e-13)
    assert np.abs(np.asarray(a) - np.asarray(b)).max() < 1e-9


def test_batch_of_detunings_equals_separate_solves():
    sd, grid = CUSTOM
    ops = build_operators(System(**copy.deepcopy(sd)))
    tlist = np.linspace(grid[0], grid[1], 61)
    rabi = lambda t: jnp.array([0.7, 0.4]) * jnp.exp(-((t - 1.2) ** 2))  # noqa: E731
    dets = np.array([[0.0, 0.0], [0.5, -0.3], [-0.8, 0.6]])
    batch = np.asarray(solve_slice(ops, rabi, dets, tlist, rtol=1e-11, atol=1e-13))
    assert batch.shape == (61, 3, 4, 4)
    for v in range(3):
        one = np.asarray(solve_slice(ops, rabi, dets[v], tlist, rtol=1e-11, atol=1e-13))
        assert np.abs(batch[:, v] - one).max() < 1e-7


def test_decay_scale_equals_scaling_the_rates():
    sd, grid = CUSTOM
    scaled = copy.deepcopy(sd)
    for d in scaled["decays"]:
        d["rate"] *= 2.5
    tlist = np.linspace(grid[0], grid[1], 61)
    rabi = lambda t: jnp.array([0.7, 0.4]) * jnp.exp(-((t - 1.2) ** 2))  # noqa: E731
    ops = build_operators(System(**copy.deepcopy(sd)))
    ops2 = build_operators(System(**scaled))
    a = solve_slice(
        ops, rabi, ops.detunings, tlist, decay_scale=2.5, rtol=1e-11, atol=1e-13
    )
    b = solve_slice(ops2, rabi, ops2.detunings, tlist, rtol=1e-11, atol=1e-13)
    assert np.abs(np.asarray(a) - np.asarray(b)).max() < 1e-9


def test_observe_stores_only_the_reduction():
    sd, grid = CUSTOM
    ops = build_operators(System(**copy.deepcopy(sd)))
    tlist = np.linspace(grid[0], grid[1], 31)
    rabi = lambda t: jnp.array([0.7, 0.4]) * jnp.exp(-((t - 1.2) ** 2))  # noqa: E731
    full = solve_slice(ops, rabi, ops.detunings, tlist, rtol=1e-11, atol=1e-13)
    pop = solve_slice(
        ops,
        rabi,
        ops.detunings,
        tlist,
        observe=lambda re, im: jnp.einsum("vii->v", re),
        rtol=1e-11,
        atol=1e-13,
    )
    assert pop.shape == (31, 1)
    assert np.allclose(
        np.asarray(pop)[:, 0], np.real(np.trace(np.asarray(full), axis1=1, axis2=2))
    )
    assert np.allclose(np.asarray(pop), 1.0, atol=1e-8)  # trace preserved


def test_no_complex_dtype_warning_from_diffrax():
    """The state is integrated as a real pair, so diffrax's complex-dtype warning
    must not fire (it flags complex support as unreliable)."""
    import warnings

    sd, grid = CUSTOM
    ops = build_operators(System(**copy.deepcopy(sd)))
    rabi = lambda t: jnp.array([0.7, 0.4]) * jnp.exp(-((t - 1.2) ** 2))  # noqa: E731
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        solve_slice(ops, rabi, ops.detunings, np.linspace(grid[0], grid[1], 21))
    assert not [w for w in caught if "omplex" in str(w.message)], [
        str(w.message) for w in caught
    ]


def test_dense_complex_dissipator_matches_matrix_exponential():
    """Non-single-element, complex collapse operators take the dense path; compare it
    with an independent matrix exponential of the Liouvillian (constant drive)."""
    from scipy.linalg import expm

    sd, _ = FIXTURES["mb_solve_lamda.json"]
    base = build_operators(System(**copy.deepcopy(sd)))
    n = base.num_states
    rng = np.random.default_rng(7)
    c_ops = 0.4 * (
        rng.normal(size=(2, n, n)) + 1j * rng.normal(size=(2, n, n))
    )  # dense and complex
    ops = dataclasses.replace(base, c_ops=c_ops, sparse=None)
    rabi_const, det = np.array([0.35, 0.2]), np.array([0.1, -0.15])
    times = np.array([0.0, 0.4, 1.1, 1.7])

    h = ops.hamiltonian(rabi_const, det)
    eye = np.eye(n)
    liou = -1j * (np.kron(eye, h) - np.kron(h.T, eye))  # column-stacking vec
    for c in c_ops:
        cdc = c.conj().T @ c
        liou += (
            np.kron(c.conj(), c) - 0.5 * np.kron(eye, cdc) - 0.5 * np.kron(cdc.T, eye)
        )
    rho0 = ops.rho0
    expected = np.array(
        [
            (expm(liou * t) @ rho0.reshape(-1, order="F")).reshape(n, n, order="F")
            for t in times
        ]
    )

    got = np.asarray(
        solve_slice(
            ops, lambda t: jnp.asarray(rabi_const), det, times, rtol=1e-12, atol=1e-14
        )
    )
    assert np.abs(got - expected).max() < 1e-9, np.abs(got - expected).max()
    assert np.abs(got[-1] - rho0).max() > 1e-2  # non-vacuous


# ---- guards ---------------------------------------------------------------------


def test_commensurate_fixed_step_is_rejected():
    sd, grid = CUSTOM
    ops = build_operators(System(**copy.deepcopy(sd)))
    tlist = np.linspace(grid[0], grid[1], 41)  # 40 equal intervals
    rabi = lambda t: jnp.zeros(2)  # noqa: E731
    with pytest.raises(ValueError, match="commensurate"):
        solve_slice(ops, rabi, ops.detunings, tlist, n_steps=2000)  # 50 per interval


def test_requires_float64():
    sd, grid = CUSTOM
    ops = build_operators(System(**copy.deepcopy(sd)))
    original = jax.config.jax_enable_x64
    try:
        jax.config.update("jax_enable_x64", False)
        with pytest.raises(RuntimeError, match="JAX_ENABLE_X64"):
            solve_slice(
                ops, lambda t: jnp.zeros(2), ops.detunings, np.linspace(0, 1, 5)
            )
    finally:
        jax.config.update("jax_enable_x64", original)


# ---- gradients --------------------------------------------------------------------

N_SAMPLES = 61
N_OUT = 41  # 40 output intervals; n_steps below is not a multiple of 40
N_FIXED = 2003


def _problem(name):
    sd, grid = ALL[name]
    ops = build_operators(System(**copy.deepcopy(sd)))
    n_f, n = len(sd["fields"]), ops.num_states
    t_samples = np.linspace(grid[0], grid[1], N_SAMPLES)
    t_out = np.linspace(grid[0], grid[1], N_OUT)
    drive = CubicDrive(t_samples)
    rng = np.random.default_rng(1)
    w = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    w = jnp.asarray(w + w.conj().T)
    x0 = np.concatenate(
        [_samples(n_f, t_samples, seed=3).ravel(), ops.detunings, [1.3]]
    )

    def loss(x, **kw):
        samples = x[: n_f * N_SAMPLES].reshape(n_f, N_SAMPLES)
        det = x[n_f * N_SAMPLES : n_f * N_SAMPLES + n_f]
        scale = x[-1]
        rho = solve_slice(ops, drive(samples), det, t_out, decay_scale=scale, **kw)
        return jnp.sum(jnp.real(jnp.einsum("ij,tji->t", w, rho)))

    # indices: field edges and interior samples, every detuning, the decay scale
    idx = []
    for k in range(n_f):
        base = k * N_SAMPLES
        idx += [
            base + i for i in (0, 1, 2, 30, N_SAMPLES - 3, N_SAMPLES - 2, N_SAMPLES - 1)
        ]
    idx += [n_f * N_SAMPLES + k for k in range(n_f)] + [n_f * N_SAMPLES + n_f]
    return loss, jnp.asarray(x0), idx


def _central_fd(f, x, idx, h):
    x = np.asarray(x, dtype=float)
    out = np.empty(len(idx))
    for j, i in enumerate(idx):
        e = np.zeros_like(x)
        e[i] = h
        out[j] = (float(f(jnp.asarray(x + e))) - float(f(jnp.asarray(x - e)))) / (2 * h)
    return out


def _rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


@pytest.mark.parametrize("name", ["custom", "mb_solve_lamda.json"])
def test_gradients_match_finite_differences(name):
    """Reverse-mode AD vs central FD on a fixed-step solve (a smooth map), for the
    field samples of every field, the detunings and the decay scale."""
    loss, x0, idx = _problem(name)
    fixed = lambda x: loss(x, n_steps=N_FIXED)  # noqa: E731
    g_fixed = np.asarray(jax.jit(jax.grad(fixed))(x0))
    fd = _central_fd(jax.jit(fixed), x0, idx, 1e-4)
    assert _rel(g_fixed[idx], fd) < 1e-5, _rel(g_fixed[idx], fd)
    assert np.abs(g_fixed[idx]).max() > 1e-4  # non-vacuous: gradients are not ~0

    # production path: adaptive solver, tight tolerance, same gradient
    adaptive = lambda x: loss(x, rtol=1e-11, atol=1e-13)  # noqa: E731
    g_pid = np.asarray(jax.jit(jax.grad(adaptive))(x0))
    assert _rel(g_pid, g_fixed) < 1e-4, _rel(g_pid, g_fixed)


def test_reverse_mode_matches_forward_mode():
    loss, x0, idx = _problem("custom")
    adaptive = lambda x: loss(x, rtol=1e-11, atol=1e-13)  # noqa: E731
    forward = lambda x: loss(  # noqa: E731
        x, rtol=1e-11, atol=1e-13, adjoint=dfx.ForwardMode()
    )
    g_rev = np.asarray(jax.jit(jax.grad(adaptive))(x0))
    g_fwd = np.asarray(jax.jit(jax.jacfwd(forward))(x0))
    assert _rel(g_rev, g_fwd) < 1e-8, _rel(g_rev, g_fwd)
