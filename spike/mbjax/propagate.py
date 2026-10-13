"""Phase 3: spatial propagation, a ``lax.scan`` over z of (Lindblad solve -> source -> AB2).

Reproduces clerq's ``PropagationSolver.solve(step="ab")`` for ``z_steps_inner = 1``:

* Omega is carried in angular units (complex); only Re(Omega) drives H (clerq's convention).
* slice 0 is driven by the *analytic* input pulse (clerq builds H_Omega from the t_func
  there); slices j >= 1 are driven by the cubic not-a-knot spline of Omega_j.
* source f_j = 2i * N(z_{j+1}) * (2 pi g) * sum_k factor_k^2 rho_{a_k b_k}, averaged over
  velocity classes with clerq's Maxwell-Boltzmann weights.
* Omega_{j+1} = Omega_j + h (1.5 f_j - 0.5 f_{j-1}), with N(z_{j+1}) used in both terms and a
  single Euler step (f_{-1} := f_0) for j = 0.

Velocity classes are carried as a leading batch axis of rho (one shared step-size
controller, no vmapped while-loops).
"""

from __future__ import annotations

import diffrax as dfx
import jax
import jax.numpy as jnp
import numpy as np

from mbjax.drive import CubicDrive
from mbjax.lindblad import build_model, check_fixed_step

jax.config.update("jax_enable_x64", True)

TWO_PI = 2.0 * np.pi


def velocity_classes(vc: dict | None, convention: str = "clerq"):
    """Detunings (cyclic, as they enter H via -2*pi*Delta*N) and normalised weights.

    convention="clerq": replicate clerq exactly. ``thermal_delta_list`` is built in
    *angular* units (2*pi*x) and then passed to ``set_H_Delta``, which multiplies by 2*pi
    again, so the Hamiltonian sees a cyclic detuning of 2*pi*x. (Suspected inconsistency in
    clerq; see PLAN.md.)  convention="physical": cyclic detuning x, i.e. thermal_delta/(2 pi).
    """
    vc = {
        "thermal_delta_min": 0.0, "thermal_delta_max": 0.0, "thermal_delta_steps": 0,
        "thermal_delta_inner_min": 0.0, "thermal_delta_inner_max": 0.0,
        "thermal_delta_inner_steps": 0, "thermal_width": 1.0, **(vc or {}),
    }
    outer = np.linspace(TWO_PI * vc["thermal_delta_min"], TWO_PI * vc["thermal_delta_max"],
                        vc["thermal_delta_steps"] + 1)
    inner = np.linspace(TWO_PI * vc["thermal_delta_inner_min"], TWO_PI * vc["thermal_delta_inner_max"],
                        vc["thermal_delta_inner_steps"] + 1)
    delta_ang = np.unique(np.concatenate([outer, inner]))
    fwhm = TWO_PI * vc["thermal_width"]
    w = 1.0 / (fwhm * np.sqrt(np.pi)) * np.exp(-((delta_ang / fwhm) ** 2))
    w = w / w.sum()
    eff = delta_ang if convention == "clerq" else delta_ang / TWO_PI
    return jnp.asarray(eff), jnp.asarray(w)


def build_propagation_model(spec: dict) -> dict:
    """Lindblad model plus the source-term index/weight arrays (field 0)."""
    m = build_model(spec)
    f = spec["fields"][0]
    facs = f.get("factors") or [1.0] * len(f["coupled_levels"])
    rows = np.array([a for a, _ in f["coupled_levels"]])
    cols = np.array([b for _, b in f["coupled_levels"]])
    m["coh_rows"], m["coh_cols"] = jnp.asarray(rows), jnp.asarray(cols)
    m["coh_w2"] = jnp.asarray(np.array(facs, dtype=float) ** 2)
    m["detuning"] = float(f.get("detuning", 0.0))
    # Sparse dissipator: if every collapse operator is a single element c_k = f_k |a_k><b_k|,
    # then c rho c^+ = |f|^2 rho_bb |a><a| and sum_k c^+c is diagonal, so the dissipator is
    # O(n^2) instead of O(k n^3). (All systems in this spike have this form.)
    cs = np.asarray(m["c"])
    nz = [np.argwhere(np.abs(c) > 0) for c in cs]
    if len(cs) and all(len(z) == 1 for z in nz):
        a = np.array([z[0][0] for z in nz]); b = np.array([z[0][1] for z in nz])
        w = np.array([abs(c[tuple(z[0])]) ** 2 for c, z in zip(cs, nz)])
        m["diss_a"], m["diss_b"], m["diss_w"] = jnp.asarray(a), jnp.asarray(b), jnp.asarray(w)
        m["diss_r"] = jnp.asarray(np.bincount(b, weights=w, minlength=m["n"]))
    return m


def _rhs_batched(t, rho, args):
    """rho: (n_v, n, n). H_v = pi*Omega(t)*X - 2*pi*(Delta_v + Delta_0)*N.

    Dissipator: dense  sum_k c_k rho c_k^+ - 1/2 {sum c^+c, rho}  (O(k n^3)), or, when every
    collapse operator is a single element f|a><b| (``diss_*`` in the model), the equivalent
    O(n^2) form: gain |a><a| (|f|^2 rho_bb) and damping -1/2 (r_i + r_j) rho_ij.
    """
    model, drive, deltas, c, cdc, gs = args
    X, N = model["X"][0], model["N"][0]
    H = jnp.pi * drive(t) * X[None] - TWO_PI * deltas[:, None, None] * N[None]
    d = -1j * (H @ rho - rho @ H)
    if "diss_a" in model:
        n = rho.shape[-1]
        pop = rho[:, model["diss_b"], model["diss_b"]]  # (n_v, k)
        gain = jnp.zeros(rho.shape[:2], rho.dtype).at[:, model["diss_a"]].add(
            gs * model["diss_w"] * pop)
        d = d + gain[:, :, None] * jnp.eye(n, dtype=rho.dtype)
        r = gs * model["diss_r"]
        return d - 0.5 * (r[:, None] + r[None, :]) * rho
    d = d + jnp.einsum("kij,vjl,kml->vim", c, rho, c.conj())
    return d - 0.5 * (cdc @ rho + rho @ cdc)


def slice_coherence(model, drive, deltas, weights, ts, *, gamma_scale=1.0, rtol=1e-9,
                    atol=1e-11, n_steps=None, adjoint=None, max_steps=1_000_000, solver=None):
    """Velocity-averaged source coherence sum_k factor_k^2 rho_{a_k b_k}(t), shape (len(ts),)."""
    check_fixed_step(ts, n_steps)
    c = model["c"] * jnp.sqrt(gamma_scale)
    cdc = model["cdc"] * gamma_scale
    n_v = deltas.shape[0]
    rho0 = jnp.broadcast_to(model["rho0"], (n_v,) + model["rho0"].shape)
    rows, cols, w2 = model["coh_rows"], model["coh_cols"], model["coh_w2"]

    def save(t, y, a):
        coh = jnp.einsum("vk,k->v", y[:, rows, cols], w2)
        return jnp.einsum("v,v->", weights.astype(coh.dtype), coh)

    sol = dfx.diffeqsolve(
        dfx.ODETerm(_rhs_batched), solver or dfx.Tsit5(), ts[0], ts[-1],
        (ts[-1] - ts[0]) / (n_steps or 1000), rho0,
        args=(model, drive, deltas + model["detuning"], c, cdc, gamma_scale),
        saveat=dfx.SaveAt(ts=ts, fn=save),
        stepsize_controller=(dfx.ConstantStepSize() if n_steps
                             else dfx.PIDController(rtol=rtol, atol=atol)),
        adjoint=adjoint or dfx.RecursiveCheckpointAdjoint(),
        max_steps=(n_steps + 2) if n_steps else max_steps,
    )
    return sol.ys


def propagate(
    model,
    omega0_ang,
    ts,
    N_z,
    g,
    h,
    deltas,
    weights,
    *,
    drive0=None,
    step="ab",
    strategy="s3",
    gamma_scale=1.0,
    rtol=1e-9,
    atol=1e-11,
    n_steps=None,
    adjoint=None,
    chunk=None,
    max_steps=1_000_000,
    solver=None,
):
    """Return Omega_zt, shape (z_steps + 1, len(ts)), angular Rabi frequency (complex).

    omega0_ang : (len(ts),) input field samples (angular units).
    drive0     : optional callable t -> cyclic Rabi frequency for slice 0 (clerq uses the
                 analytic pulse there); default is the spline of omega0_ang.
    N_z        : (z_steps,) number density at z_{j+1}.   g: cyclic interaction strength.
    """
    spline = CubicDrive(ts)
    kw = dict(gamma_scale=gamma_scale, rtol=rtol, atol=atol, n_steps=n_steps, adjoint=adjoint,
              max_steps=max_steps, solver=solver)
    g_ang = TWO_PI * g

    def coh_of(drive):
        return slice_coherence(model, drive, deltas, weights, ts, **kw)

    def advance(omega, coh, coh_prev, N):
        if step == "euler":
            return omega + h * 2j * N * g_ang * coh
        return omega + h * 2j * N * g_ang * (1.5 * coh - 0.5 * coh_prev)

    omega0 = jnp.asarray(omega0_ang).astype(jnp.complex128)
    # slice 0 (peeled): analytic drive if given, Euler bootstrap
    d0 = drive0 if drive0 is not None else spline(jnp.real(omega0) / TWO_PI)
    coh0 = coh_of(lambda t: d0(t))
    omega1 = omega0 + h * 2j * N_z[0] * g_ang * coh0

    def body(carry, N):
        omega, coh_prev = carry
        drive = spline(jnp.real(omega) / TWO_PI)
        coh = coh_of(drive)
        omega_next = advance(omega, coh, coh_prev, N)
        return (omega_next, coh), omega_next

    if strategy in ("s3", "s5"):
        body = jax.checkpoint(body)

    xs = N_z[1:]
    if strategy == "s5" and xs.shape[0] > 0:
        L = xs.shape[0]
        c = chunk or max(1, int(round(np.sqrt(L))))
        n_chunks, rem = divmod(L, c)

        @jax.checkpoint
        def chunk_body(carry, xs_c):
            return jax.lax.scan(body, carry, xs_c)

        carry, outs = jax.lax.scan(chunk_body, (omega1, coh0), xs[: n_chunks * c].reshape(n_chunks, c))
        outs = outs.reshape(n_chunks * c, -1)
        if rem:
            carry, outs_r = jax.lax.scan(body, carry, xs[n_chunks * c:])
            outs = jnp.concatenate([outs, outs_r])
    else:
        _, outs = jax.lax.scan(body, (omega1, coh0), xs)
    return jnp.concatenate([omega0[None], omega1[None], outs])
