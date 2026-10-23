"""Operator data of a :class:`~clerq.system.System` for the autodiff solvers.

Everything here is plain NumPy: the arrays are extracted once from a ``System``
and converted to JAX arrays by the solver. Conventions are those of ``System``
(frequencies are cyclic; the Hamiltonian carries ``2 pi`` and ``pi`` factors)::

    H(t) = H_0 + sum_k [ pi * Omega_k(t) * X_k  -  s_k * 2 pi * Delta_k * N_k ]

where field ``k`` couples the pairs ``(a, b)`` with signed factors ``f``:
``X_k = sum f (|a><b| + |b><a|)``, ``N_k`` projects on the upper levels the field
couples, and ``s_k = +1`` (``-1`` if ``detuning_positive`` is False). The source
term of the Maxwell equation weights the coherence of each pair by its signed
factor: ``sum f rho_{ab}``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from clerq.system import System

_TOL = 1e-14


@dataclass(frozen=True)
class SparseDissipator:
    """Dissipator of collapse operators that are all single elements ``f |a><b|``.

    Then ``c rho c^+ = |f|^2 rho_bb |a><a|`` and ``sum c^+ c`` is diagonal, so the
    dissipator costs ``O(n^2)`` rather than ``O(k n^3)``::

        D(rho)_ij = delta_ij sum_{k: a_k = i} w_k rho_{b_k b_k}
                    - (r_i + r_j) / 2 * rho_ij

    Attributes:
        a: target level of each operator, shape ``(k,)``.
        b: source level of each operator, shape ``(k,)``.
        w: ``|f|^2`` of each operator (rates included), shape ``(k,)``.
        r: total rate out of each level, ``r_i = sum_{k: b_k = i} w_k``, ``(n,)``.
    """

    a: np.ndarray
    b: np.ndarray
    w: np.ndarray
    r: np.ndarray

    def apply(self, rho: np.ndarray) -> np.ndarray:
        """Return ``D(rho)`` for ``rho`` of shape ``(..., n, n)`` (NumPy reference)."""
        rho = np.asarray(rho)
        n = rho.shape[-1]
        out = -0.5 * (self.r[:, None] + self.r[None, :]) * rho
        gain = np.zeros(rho.shape[:-2] + (n,), dtype=rho.dtype)
        np.add.at(gain, (Ellipsis, self.a), self.w * rho[..., self.b, self.b])
        idx = np.arange(n)
        out[..., idx, idx] += gain
        return out


@dataclass(frozen=True)
class DecayChannels:
    """The decay channels of a ``System``, one entry per ``(decay, channel)`` pair.

    Every channel is a collapse operator ``f sqrt(2 pi rate) |a><b|``, so it is a
    single element and the dissipator is linear in its weight
    ``w_k = 2 pi rate f^2``. Keeping the channels (rather than the finished
    operators) makes each decay's ``rate`` a parameter, differentiable at any
    value including zero.

    Attributes:
        a: target level of each channel, shape ``(k,)``.
        b: source level of each channel, shape ``(k,)``.
        f2: squared factor ``f^2`` of each channel, shape ``(k,)``.
        idx: index of the decay (entry of ``system.decays``) each channel belongs to.
        rates: cyclic ``rate`` of each decay, shape ``(n_decays,)``.
    """

    a: np.ndarray
    b: np.ndarray
    f2: np.ndarray
    idx: np.ndarray
    rates: np.ndarray


@dataclass(frozen=True)
class SystemOperators:
    """Operators and field data extracted from a ``System``.

    Attributes:
        num_states: dimension ``n`` of the system.
        H0: bare Hamiltonian ``2 pi diag(energies)``, shape ``(n, n)``, real.
        X: coupling operators, shape ``(n_fields, n, n)``, real symmetric, so that
            the field's interaction is ``pi * Omega_k * X_k`` (``Omega`` cyclic).
        Y: antisymmetric coupling operators, shape ``(n_fields, n, n)``, so that a
            complex Rabi frequency ``Omega`` drives ``pi (Re Omega X_k + i Im Omega Y_k)``
            (``Omega`` on ``rho_ab``'s Hamiltonian element, ``conj(Omega)`` on the
            transpose). Built from ``pairs`` when not given.
        N: detuning operators, shape ``(n_fields, n, n)``, real diagonal.
        detuning_sign: ``s_k`` of each field (``+1`` or ``-1``).
        detunings: fixed cyclic detuning of each field, shape ``(n_fields,)``.
        doppler_factors: multiplier of the thermal detuning for each field.
        rabi_freqs: peak cyclic Rabi frequency of each field.
        pairs: per field ``(rows, cols, factors)`` of the coupled coherences.
        c_ops: nonzero collapse operators, shape ``(k, n, n)`` (``k`` may be 0).
        sparse: the O(n^2) form of the dissipator if every collapse operator is a
            single element, else ``None``.
        rho0: initial density matrix, shape ``(n, n)``.
        Z: Zeeman operator ``diag(g_F m_F)``, shape ``(n, n)``; the Zeeman part
            of ``H0`` is ``2 pi b_field Z``.
        b_field: the system's cyclic B-field in natural units (already in ``H0``).
        has_b_field: whether the system has a ``b_field`` block (``Z`` and ``b_field``
            default to no Zeeman term).
    """

    num_states: int
    H0: np.ndarray
    X: np.ndarray
    N: np.ndarray
    detuning_sign: np.ndarray
    detunings: np.ndarray
    doppler_factors: np.ndarray
    rabi_freqs: np.ndarray
    pairs: tuple[tuple[np.ndarray, np.ndarray, np.ndarray], ...]
    c_ops: np.ndarray
    sparse: SparseDissipator | None
    rho0: np.ndarray
    Z: np.ndarray | None = None
    b_field: float = 0.0
    has_b_field: bool = False
    Y: np.ndarray | None = None
    decays: DecayChannels | None = None

    def __post_init__(self):
        if self.Y is None:
            Y = np.zeros_like(self.X)
            for k, (rows, cols, factors) in enumerate(self.pairs):
                np.add.at(Y[k], (rows, cols), factors)
                np.add.at(Y[k], (cols, rows), -factors)
            object.__setattr__(self, "Y", Y)
        if self.Z is None:
            object.__setattr__(self, "Z", np.zeros((self.num_states, self.num_states)))

    @property
    def num_fields(self) -> int:
        return self.X.shape[0]

    def hamiltonian(
        self, rabi_freqs: np.ndarray, detunings: np.ndarray | None = None
    ) -> np.ndarray:
        """``H`` for constant cyclic (possibly complex) Rabi frequencies and detunings.

        NumPy reference.
        """
        det = self.detunings if detunings is None else np.asarray(detunings)
        rabi = np.asarray(rabi_freqs)
        return (
            self.H0
            + np.pi * np.einsum("k,kij->ij", rabi.real, self.X)
            + 1j * np.pi * np.einsum("k,kij->ij", rabi.imag, self.Y)
            - 2.0 * np.pi * np.einsum("k,kij->ij", self.detuning_sign * det, self.N)
        )

    def source_coherence(self, rho: np.ndarray) -> np.ndarray:
        """``sum_pairs factor * rho_ab`` per field, shape ``(n_fields, ...)``.

        ``rho`` has shape ``(..., n, n)``; this is the quantity the Maxwell source
        term is proportional to (NumPy reference).
        """
        rho = np.asarray(rho)
        out = np.zeros((self.num_fields,) + rho.shape[:-2], dtype=complex)
        for k, (rows, cols, factors) in enumerate(self.pairs):
            out[k] = np.einsum("...p,p->...", rho[..., rows, cols], factors)
        return out


def _sparse_dissipator(c_ops: np.ndarray, n: int) -> SparseDissipator | None:
    """Return the sparse form if every collapse operator is a single element."""
    if len(c_ops) == 0:
        return None
    a, b, w = [], [], []
    for c in c_ops:
        nz = np.argwhere(np.abs(c) > _TOL)
        if len(nz) != 1:
            return None
        i, j = nz[0]
        a.append(i)
        b.append(j)
        w.append(abs(c[i, j]) ** 2)
    a, b, w = np.array(a), np.array(b), np.array(w)
    return SparseDissipator(a=a, b=b, w=w, r=np.bincount(b, weights=w, minlength=n))


def _decay_channels(system: System) -> DecayChannels | None:
    """The channels of ``system.decays`` (``None`` if it has no decays)."""
    a, b, f2, idx = [], [], [], []
    for i, d in enumerate(system.decays):
        channels = d["channels"]
        factors = d.get("factors", [1.0] * len(channels))
        for (to, frm), f in zip(channels, factors, strict=True):
            a.append(to)
            b.append(frm)
            f2.append(float(f) ** 2)
            idx.append(i)
    if not idx:
        return None
    return DecayChannels(
        a=np.array(a, dtype=int),
        b=np.array(b, dtype=int),
        f2=np.array(f2),
        idx=np.array(idx, dtype=int),
        rates=np.array([float(d["rate"]) for d in system.decays]),
    )


def build_operators(system: System) -> SystemOperators:
    """Extract the operators of a ``System`` for the autodiff solvers.

    Args:
        system: the system, with its fields and decays.

    Returns:
        The operator data.

    Raises:
        NotImplementedError: if a field is counter-propagating.
    """
    n = system.num_states
    for f in system.fields:
        if f.counter_propagating:
            raise NotImplementedError(
                f"Field {f.label!r} is counter-propagating; clerq.autodiff supports "
                "co-propagating fields only."
            )

    H0 = np.asarray(system.H_0.full())
    if np.abs(H0.imag).max(initial=0.0) > _TOL:
        raise ValueError("The bare Hamiltonian H_0 must be real.")
    H0 = H0.real

    n_f = len(system.fields)
    X = np.zeros((n_f, n, n))
    N = np.zeros((n_f, n, n))
    pairs = []
    for k, f in enumerate(system.fields):
        rows = np.array([c[0] for c in f.coupled_levels], dtype=int)
        cols = np.array([c[1] for c in f.coupled_levels], dtype=int)
        factors = np.asarray(f.factors, dtype=float)
        np.add.at(X[k], (rows, cols), factors)
        np.add.at(X[k], (cols, rows), factors)
        for i in f.upper_levels():
            N[k, i, i] = 1.0
        pairs.append((rows, cols, factors))

    c_ops = np.array([np.asarray(c.full()) for c in system.c_ops]).reshape(-1, n, n)
    c_ops = c_ops[np.abs(c_ops).reshape(len(c_ops), -1).max(axis=1) > 0]

    return SystemOperators(
        num_states=n,
        H0=H0,
        X=X,
        N=N,
        detuning_sign=np.array(
            [1.0 if f.detuning_positive else -1.0 for f in system.fields]
        ),
        detunings=np.array([f.detuning for f in system.fields], dtype=float),
        doppler_factors=np.array([f.factor_doppler_shift for f in system.fields]),
        rabi_freqs=np.array([f.rabi_freq for f in system.fields], dtype=float),
        pairs=tuple(pairs),
        c_ops=c_ops,
        sparse=_sparse_dissipator(c_ops, n),
        rho0=np.asarray(system.initial_state.full()),
        Z=np.diag(system.zeeman_coeffs),
        b_field=float(system.b_field_natural),
        has_b_field=system.b_field is not None,
        decays=_decay_channels(system),
    )
