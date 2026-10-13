"""T0.3: system specs are valid, physical, and agree with clerq's own operators (clerq.System)."""

import copy

import numpy as np
import pytest
from clerq import System

from mbjax.systems import DIMS, liouvillian, make_system, operators

NAMES = list(DIMS)


@pytest.mark.parametrize("name", NAMES)
def test_dimension(name):
    assert make_system(name)["num_states"] == DIMS[name]


@pytest.mark.parametrize("name", NAMES)
def test_branching_sums_to_one(name):
    spec = make_system(name)
    d = spec["decays"][0]
    n_g = len({a for a, _ in d["channels"]})
    excited = sorted({b for _, b in d["channels"]})
    for e in excited:
        tot = sum(f**2 for (a, b), f in zip(d["channels"], d["factors"]) if b == e)
        assert tot == pytest.approx(1.0, abs=1e-12), (name, e, n_g)


@pytest.mark.parametrize("name", NAMES)
def test_liouvillian_trace_preserving_and_physical(name):
    ops = operators(make_system(name))
    L = liouvillian(ops, rabi=[0.7], detuning=[0.3])
    n = ops["n"]
    tr = np.eye(n).reshape(-1, order="F")  # <vec(I)| picks the trace
    assert np.allclose(tr @ L, 0.0, atol=1e-10)  # d/dt Tr rho = 0
    ev = np.linalg.eigvals(L)
    assert ev.real.max() < 1e-9  # no growing modes
    assert np.sum(np.abs(ev) < 1e-9) >= 1  # a steady state exists


@pytest.mark.parametrize("name", NAMES)
def test_matches_clerq_operators(name):
    spec = make_system(name)
    ops = operators(spec)
    atom = System(**copy.deepcopy(spec))
    # collapse operators: same count, same matrices, same order
    assert len(atom.c_ops) == len(ops["c_ops"])
    for qc, c in zip(atom.c_ops, ops["c_ops"], strict=True):
        assert np.allclose(qc.full(), c, atol=1e-12)
    # H_Omega (field is time-dependent only if a t_func is set; here it is constant)
    for k, f in enumerate(atom.fields):
        Hk = atom.H_Omega_list[k]
        Hk = Hk if not isinstance(Hk, list) else Hk[0]
        assert np.allclose(Hk.full(), np.pi * f.rabi_freq * ops["X"][k], atol=1e-12)
    # H_Delta at detuning 0.37
    atom.set_H_Delta([0.37] * len(atom.fields))
    ref = -2 * np.pi * 0.37 * sum(ops["N"])
    assert np.allclose(atom.H_Delta.full(), ref, atol=1e-12)
    # initial state
    assert np.allclose(atom.initial_state.full(), ops["rho0"], atol=1e-12)
