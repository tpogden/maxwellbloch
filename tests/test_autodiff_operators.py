"""The operators extracted for clerq.autodiff equal those of ``System`` itself.

``System`` is the live oracle: every comparison is against the object clerq builds
from the same problem definition, so a change to a clerq convention fails here.
"""

import copy
import glob
import json
import os

import numpy as np
import pytest

pytest.importorskip("jax")
pytest.importorskip("diffrax")

from clerq.autodiff._operators import (  # noqa: E402
    _sparse_dissipator,
    build_operators,
)
from clerq.system import System  # noqa: E402

pytestmark = pytest.mark.jax

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")


def _system_dict(path):
    with open(path) as f:
        d = json.load(f)
    if "system" in d:
        return d["system"]
    if "num_states" in d:
        return d
    return None


FIXTURES = {
    os.path.basename(p): sd
    for p in sorted(glob.glob(os.path.join(JSON_DIR, "*.json")))
    if (sd := _system_dict(p)) is not None
}

CUSTOM = {
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
        },
        {
            "coupled_levels": [[2, 3]],
            "detuning": -0.2,
            "factor_doppler_shift": 0.5,
            "rabi_freq": 0.9,
        },
    ],
    "initial_state": [0.5, 0.5, 0.0, 0.0],
}
ALL = {**FIXTURES, "custom": CUSTOM}


def _dense_dissipator(c_ops, rho):
    """sum_k c rho c^+ - 1/2 {sum_k c^+ c, rho} for rho of shape (..., n, n)."""
    out = np.zeros_like(rho)
    cdc = np.zeros(rho.shape[-2:], dtype=complex)
    for c in c_ops:
        out += c @ rho @ c.conj().T
        cdc += c.conj().T @ c
    return out - 0.5 * (cdc @ rho + rho @ cdc)


def _random_states(n, count, seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(count, n, n)) + 1j * rng.normal(size=(count, n, n))
    return a @ a.conj().transpose(0, 2, 1)


def test_oracle_set_is_not_vacuous():
    """The comparisons below must actually cover the structure they claim to."""
    systems = [System(**copy.deepcopy(d)) for d in ALL.values()]
    assert len(FIXTURES) >= 10
    assert any(len(s.fields) >= 2 for s in systems), "no multi-field system"
    assert any(any(s.energies) for s in systems if s.energies), "no bare energies"
    assert any(len(s.decays) == 0 for s in systems), "no decay-free system"
    assert any(f != 1.0 for s in systems for fl in s.fields for f in fl.factors)
    assert any(not fl.detuning_positive for s in systems for fl in s.fields)


@pytest.mark.parametrize("name", list(ALL))
class TestMatchesSystem:
    def _build(self, name):
        system = System(**copy.deepcopy(ALL[name]))
        return system, build_operators(system)

    def test_dimension_and_initial_state(self, name):
        system, ops = self._build(name)
        assert ops.num_states == system.num_states
        assert np.allclose(ops.rho0, system.initial_state.full(), atol=1e-14)

    def test_bare_hamiltonian(self, name):
        system, ops = self._build(name)
        assert np.allclose(ops.H0, system.H_0.full(), atol=1e-12)

    def test_interaction_hamiltonian(self, name):
        system, ops = self._build(name)
        assert ops.num_fields == len(system.fields)
        for k, f in enumerate(system.fields):
            h = system.H_Omega_list[k]
            h = h[0] if isinstance(h, list) else h  # [operator, time function]
            assert np.allclose(h.full(), np.pi * f.rabi_freq * ops.X[k], atol=1e-12)

    def test_detuning_hamiltonian(self, name):
        system, ops = self._build(name)
        test = 0.37 * (1.0 + np.arange(ops.num_fields))
        h_delta = system.set_H_Delta(list(test)).full()
        ref = -2.0 * np.pi * np.einsum("k,kij->ij", ops.detuning_sign * test, ops.N)
        assert np.allclose(h_delta, ref, atol=1e-12)

    def test_field_scalars(self, name):
        system, ops = self._build(name)
        assert np.allclose(ops.detunings, [f.detuning for f in system.fields])
        assert np.allclose(
            ops.doppler_factors, [f.factor_doppler_shift for f in system.fields]
        )
        assert np.allclose(ops.rabi_freqs, [f.rabi_freq for f in system.fields])

    def test_collapse_operators(self, name):
        system, ops = self._build(name)
        nonzero = [c.full() for c in system.c_ops if np.abs(c.full()).max() > 0]
        assert len(ops.c_ops) == len(nonzero)
        for mine, theirs in zip(ops.c_ops, nonzero, strict=True):
            assert np.allclose(mine, theirs, atol=1e-14)

    def test_sparse_dissipator_equals_dense(self, name):
        system, ops = self._build(name)
        if len(ops.c_ops) == 0:
            assert ops.sparse is None
            return
        assert ops.sparse is not None  # System collapse operators are single elements
        rho = _random_states(ops.num_states, 3, seed=1)
        expected = _dense_dissipator(ops.c_ops, rho)
        scale = max(1.0, np.abs(expected).max())
        assert np.abs(ops.sparse.apply(rho) - expected).max() < 1e-12 * scale

    def test_source_coherence_weights(self, name):
        """The Maxwell source term weights each coherence by its signed factor."""
        system, ops = self._build(name)
        if not system.fields:
            pytest.skip("no fields")
        states = _random_states(ops.num_states, 4, seed=2)
        expected = system.get_fields_sum_coherence(states_t=states)
        assert np.allclose(ops.source_coherence(states), expected, atol=1e-12)

    def test_full_hamiltonian(self, name):
        """H_0 + H_Delta + H_Omega assembled by hand equals clerq's static sum."""
        system, ops = self._build(name)
        if not system.fields:
            pytest.skip("no fields")
        h_omega = sum(
            (h[0] if isinstance(h, list) else h).full() for h in system.H_Omega_list
        )
        total = system.H_0.full() + system.H_Delta.full() + h_omega
        mine = ops.hamiltonian(ops.rabi_freqs, ops.detunings)
        assert np.allclose(mine, total, atol=1e-12)


def test_dense_collapse_operator_has_no_sparse_form():
    c = np.array([[[0.0, 1.0], [1.0, 0.0]]])  # two nonzero entries
    assert _sparse_dissipator(c, 2) is None


def test_counter_propagating_field_is_rejected():
    spec = copy.deepcopy(CUSTOM)
    spec["fields"][0]["counter_propagating"] = True
    spec["fields"][0].pop("factor_doppler_shift", None)
    with pytest.raises(NotImplementedError, match="co-propagating"):
        build_operators(System(**spec))
