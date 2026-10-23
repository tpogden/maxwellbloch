"""Independent NumPy/SciPy reference for the propagation problem with an optional full complex (RWA) drive.
real_drive=True reproduces the deprecated convention (only Re(Omega) enters H); False uses Omega on one
off-diagonal and Omega* on the other, as clerq does by default."""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from clerq.autodiff._operators import build_operators


def reference(
    s,
    omega_in_cyc,
    real_drive,
    zs=30,
    zi=20,
    rtol=1e-7,
    atol=1e-10,
    detunings=None,
    verbose=False,
):
    ops = build_operators(s.system)
    n = ops.num_states
    nf = len(ops.pairs)
    tl = np.asarray(s.tlist)
    zl = np.linspace(s.z_min, s.z_max, zs + 1)
    g = (
        2 * np.pi * np.asarray(s.interaction_strengths, float)
    )  # clerq: g = 2 pi * interaction_strength
    thermal = np.asarray(s.thermal_delta_list, float) / (2 * np.pi)
    wts = np.asarray(s.thermal_weights, float)
    wts = wts / wts.sum()
    V = len(thermal)
    det = np.asarray(ops.detunings if detunings is None else detunings, float)
    shifts = (
        thermal[:, None] * np.asarray(ops.doppler_factors)[None, :]
    )  # (V, nf) cyclic
    Hs = np.stack(
        [
            ops.H0
            - 2
            * np.pi
            * np.einsum("k,kij->ij", (det + shifts[v]) * ops.detuning_sign, ops.N)
            for v in range(V)
        ]
    ).astype(complex)
    cops = np.asarray(ops.c_ops, complex)
    CdC = (
        np.einsum("kji,kjl->il", cops.conj(), cops)
        if len(cops)
        else np.zeros((n, n), complex)
    )
    rho0 = np.broadcast_to(np.asarray(ops.rho0, complex), (V, n, n)).copy()

    def dens(z):  # square density on [0, 1] (clerq default)
        return 1.0 if (0.0 <= z <= 1.0) else 0.0

    def sumcoh(Om_cyc):
        splines = [CubicSpline(tl, Om_cyc[k]) for k in range(nf)]

        def rhs(t, y):
            rho = y.reshape(V, n, n)
            H = Hs.copy()
            for k in range(nf):
                o = splines[k](t)
                o = complex(o.real) if real_drive else complex(o)
                rows, cols, facs = ops.pairs[k]
                for r, c, f in zip(rows, cols, facs, strict=True):
                    H[:, r, c] += np.pi * o * f
                    H[:, c, r] += np.pi * np.conj(o) * f
            d = -1j * (H @ rho - rho @ H)
            if len(cops):
                d += np.einsum("kil,vlm,knm->vin", cops, rho, cops.conj()) - 0.5 * (
                    CdC[None] @ rho + rho @ CdC[None]
                )
            return d.reshape(-1)

        sol = solve_ivp(
            rhs,
            (tl[0], tl[-1]),
            rho0.reshape(-1),
            t_eval=tl,
            method="RK45",
            rtol=rtol,
            atol=atol,
        )
        R = sol.y.reshape(V, n, n, -1)  # (V, n, n, nt)
        Ravg = np.tensordot(wts, R, axes=(0, 0))  # (n, n, nt)
        out = np.zeros((nf, len(tl)), complex)
        for k in range(nf):
            rows, cols, facs = ops.pairs[k]
            for r, c, f in zip(rows, cols, facs, strict=True):
                out[k] += f * Ravg[r, c]
        return out

    Om = np.asarray(omega_in_cyc, complex) * 2 * np.pi  # angular, as clerq
    h = (zl[1] - zl[0]) / zi
    z = zl[0]
    prev = None
    out = [Om.copy()]
    for j in range(zs):
        for _ in range(zi):
            znext = z + h
            coh = sumcoh(Om / (2 * np.pi))
            d_this = 2j * g[:, None] * coh
            if prev is None:
                step = d_this
            else:
                step = 1.5 * d_this - 0.5 * prev
            Om = Om + h * dens(znext) * step
            prev = d_this
            z = znext
        out.append(Om.copy())
        if verbose:
            print("slice", j + 1, "/", zs, flush=True)
    return np.array(out)  # (zs+1, nf, nt) angular
