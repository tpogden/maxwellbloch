"""T0.3: shared system specs for the spike, in clerq's atom-dict format.

A spec is a plain dict accepted by ``clerq.System(**spec)`` (and so by clerq's JSON
files): ``num_states``, ``decays``, ``fields``, ``initial_state``. ``operators`` turns
a spec into the dense operator pieces in clerq's conventions, without QuTiP:

    H(t) = -(+/-)2*pi*Delta_k * N_k  +  pi * Omega_k(t) * X_k      (per field k)
    c_i  = factor * sqrt(2*pi*rate) * |a><b|                       (per channel)

so frequencies are cyclic (rate Gamma in units where Gamma = 1 means 2*pi*Gamma rad).

Systems (one linearly x-polarised probe field coupling the whole manifold; ground
levels first in ascending m, then excited):

    two_level   dim  2   |g> <-> |e>
    j12_12      dim  4   Jg=1/2 -> Je=1/2
    j1_2_res    dim  9   Jg=1 -> Je=2 (sigma+), plus a metastable reservoir level (index 8)
                         receiving 10% of the decay from every excited level
    j72_72      dim 16   Jg=7/2 -> Je=7/2
"""

from __future__ import annotations

from math import pi, sqrt

import numpy as np
from clerq.angmom import wigner_3j


def _cg(jg, mg, je, me, q):
    """<jg mg; 1 q | je me> from the 3j symbol (zero if m's don't conserve)."""
    if mg + q != me:
        return 0.0
    return (-1) ** (jg - 1 + me) * sqrt(2 * je + 1) * wigner_3j(jg, 1, je, mg, q, -me)


def _two_level(rabi_freq: float = 1.0, gamma: float = 1.0) -> dict:
    return {
        "label": "two_level",
        "num_states": 2,
        "decays": [{"rate": gamma, "channels": [[0, 1]], "factors": [1.0]}],
        "initial_state": [1.0, 0.0],
        "fields": [
            {
                "label": "probe",
                "coupled_levels": [[0, 1]],
                "factors": [1.0],
                "detuning": 0.0,
                "rabi_freq": rabi_freq,
            }
        ],
    }


def _j_to_j(
    jg: float,
    je: float,
    label: str,
    rabi_freq=1.0,
    gamma=1.0,
    reservoir: float = 0.0,
    pol: str = "x",
) -> dict:
    """Probe on a Jg -> Je manifold: ``pol="x"`` is linear (q=+/-1 with -1/sqrt2, +1/sqrt2),
    ``pol="sigma+"`` is circular (q=+1 only, factor = CG).

    Note for the propagation tests: clerq's field source weights each coupled pair by
    factor**2 while H uses factor, so the medium's response goes as sum(factor**3). For the
    linear polarisation with a symmetric manifold that sum vanishes identically (J=1 -> 2),
    which makes a propagation test vacuous; use ``sigma+`` there.

    If ``reservoir`` > 0, one extra level (last index) is appended, uncoupled to the
    field and fed by that fraction of every excited level's decay."""
    mgs = [-jg + i for i in range(int(round(2 * jg)) + 1)]
    mes = [-je + i for i in range(int(round(2 * je)) + 1)]
    n_g = len(mgs)
    g_idx = {m: i for i, m in enumerate(mgs)}
    e_idx = {m: n_g + i for i, m in enumerate(mes)}

    coupled, factors = [], []
    for mg in mgs:
        for q, s in ((+1, -1.0), (-1, +1.0)) if pol == "x" else ((+1, 1.0),):
            me = mg + q
            if me in e_idx:
                c = _cg(jg, mg, je, me, q)
                if abs(c) > 1e-12:
                    coupled.append([g_idx[mg], e_idx[me]])
                    factors.append(s / sqrt(2) * c if pol == "x" else c)

    # spontaneous decay e -> g: one collapse operator per (g, e); branching b = CG^2
    # with q = me - mg in {-1, 0, +1}, so sum_g b = 1 for every excited level.
    channels, dfactors = [], []
    for me in mes:
        for mg in mgs:
            q = me - mg
            if abs(q) <= 1:
                b = _cg(jg, mg, je, me, q) ** 2
                if b > 1e-12:
                    channels.append([g_idx[mg], e_idx[me]])
                    dfactors.append(sqrt(b * (1.0 - reservoir)))

    n = n_g + len(mes)
    if reservoir > 0:
        for me in mes:
            channels.append([n, e_idx[me]])
            dfactors.append(sqrt(reservoir))
        n += 1
    pop = [1.0 / n_g] * n_g + [0.0] * (n - n_g)  # unpolarised ground manifold
    return {
        "label": label,
        "num_states": n,
        "decays": [{"rate": gamma, "channels": channels, "factors": dfactors}],
        "initial_state": pop,
        "fields": [
            {
                "label": "probe",
                "coupled_levels": coupled,
                "factors": factors,
                "detuning": 0.0,
                "rabi_freq": rabi_freq,
            }
        ],
    }


SYSTEMS = {
    "two_level": lambda **kw: _two_level(**kw),
    "j12_12": lambda **kw: _j_to_j(0.5, 0.5, "j12_12", **kw),
    "j1_2_res": lambda **kw: _j_to_j(1.0, 2.0, "j1_2_res", reservoir=0.1, pol="sigma+", **kw),
    "j72_72": lambda **kw: _j_to_j(3.5, 3.5, "j72_72", **kw),
    # linear polarisation on J=1 -> 2: its factors are antisymmetric, so the response vanished
    # under clerq <= 0.13.1 (factor**2 source weights); fine since the dipole-weight fix.
    "j1_2_res_x": lambda **kw: _j_to_j(1.0, 2.0, "j1_2_res_x", reservoir=0.1, pol="x", **kw),
}
DIMS = {"two_level": 2, "j12_12": 4, "j1_2_res": 9, "j72_72": 16, "j1_2_res_x": 9}


def make_system(name: str, **kw) -> dict:
    """Return a fresh spec dict (safe to mutate / pass to System)."""
    return SYSTEMS[name](**kw)


def _sigma(n: int, a: int, b: int) -> np.ndarray:
    s = np.zeros((n, n), dtype=complex)
    s[a, b] = 1.0
    return s


def operators(spec: dict) -> dict:
    """Dense operator pieces for a spec, in clerq's conventions (no QuTiP).

    Returns dict with:
        n        : num_states
        X        : [n x n] per field, so H_Omega_k = pi * Omega_k(t) * X[k]
        N        : [n x n] per field, so H_Delta_k = -sgn*2*pi*Delta_k * N[k]
        c_ops    : list of n x n collapse operators (rates included)
        rho0     : n x n initial density matrix (diagonal)
    """
    n = spec["num_states"]
    X, N = [], []
    for f in spec["fields"]:
        facs = f.get("factors") or [1.0] * len(f["coupled_levels"])
        X.append(
            sum(
                fac * (_sigma(n, a, b) + _sigma(n, b, a))
                for (a, b), fac in zip(f["coupled_levels"], facs, strict=True)
            )
        )
        N.append(
            sum(_sigma(n, i, i) for i in sorted({b for _, b in f["coupled_levels"]}))
        )
    c_ops = []
    for d in spec.get("decays", []):
        facs = d.get("factors") or [1.0] * len(d["channels"])
        for (a, b), fac in zip(d["channels"], facs, strict=True):
            c_ops.append(fac * sqrt(2 * pi * d["rate"]) * _sigma(n, a, b))
    pop = spec.get("initial_state") or [1.0] + [0.0] * (n - 1)
    return {
        "n": n,
        "X": X,
        "N": N,
        "c_ops": c_ops,
        "rho0": np.diag(np.asarray(pop, dtype=float)).astype(complex),
    }


def liouvillian(ops: dict, rabi: list[float], detuning: list[float]) -> np.ndarray:
    """Dense n^2 x n^2 Liouvillian (column-stacking vec) at constant Omega_k, Delta_k,
    for cross-checks. H = sum_k [-2*pi*Delta_k N_k + pi*Omega_k X_k]."""
    n = ops["n"]
    H = np.zeros((n, n), dtype=complex)
    for Xk, Nk, om, de in zip(ops["X"], ops["N"], rabi, detuning, strict=True):
        H += -2 * pi * de * Nk + pi * om * Xk
    I = np.eye(n)
    # vec(A rho B) = (B^T kron A) vec(rho), column-stacking
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for c in ops["c_ops"]:
        cdc = c.conj().T @ c
        L += np.kron(c.conj(), c) - 0.5 * np.kron(I, cdc) - 0.5 * np.kron(cdc.T, I)
    return L
