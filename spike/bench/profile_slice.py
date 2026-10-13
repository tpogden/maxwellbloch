"""Why is the stretch point slow? Time ONE z-slice (dim 9, n_t 2000, n_v 50) under variants."""

import json
import sys
import time

import diffrax as dfx
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

dim, n_t, n_v = 9, 2000, 50
n = dim * dim
k = jax.random.split(jax.random.PRNGKey(0), 3)
B = 0.3 * jax.random.normal(k[0], (n, n))
C = 0.3 * jax.random.normal(k[1], (n, n))
C = C - C.T - 0.5 * jnp.eye(n)
D = 0.3 * jax.random.normal(k[2], (n, n))
D = D - D.T
ts = jnp.linspace(0.0, 1.0, n_t)
deltas = jnp.linspace(-1.0, 1.0, n_v)
y0 = jnp.zeros(n).at[0].set(1.0)
field0 = jnp.exp(-((ts - 0.5) ** 2) / 0.02)


def solve(vf, y0_, saveat, adjoint=None):
    return dfx.diffeqsolve(
        dfx.ODETerm(vf), dfx.Tsit5(), 0.0, 1.0, 1.0 / (n_t - 1), y0_,
        saveat=saveat, stepsize_controller=dfx.ConstantStepSize(),
        adjoint=adjoint or dfx.RecursiveCheckpointAdjoint(), max_steps=n_t + 8,
    )


def base(field):  # vmap over classes, save full state
    def cls(delta):
        vf = lambda t, y, a: (C + jnp.interp(t, ts, field) * B + delta * D) @ y
        return solve(vf, y0, dfx.SaveAt(ts=ts)).ys[:, 0]
    return jax.vmap(cls)(deltas).mean(0)


def pfn(field):  # vmap over classes, save only polarisation
    def cls(delta):
        vf = lambda t, y, a: (C + jnp.interp(t, ts, field) * B + delta * D) @ y
        sa = dfx.SaveAt(ts=ts, fn=lambda t, y, a: y[0])
        return solve(vf, y0, sa).ys
    return jax.vmap(cls)(deltas).mean(0)


def mat(field):  # no vmap: state is (n, n_v), three GEMMs, save only polarisation
    Y0 = jnp.tile(y0[:, None], (1, n_v))
    def vf(t, Y, a):
        E = jnp.interp(t, ts, field)
        return C @ Y + E * (B @ Y) + (D @ Y) * deltas[None, :]
    sa = dfx.SaveAt(ts=ts, fn=lambda t, Y, a: Y[0])
    return solve(vf, Y0, sa).ys.mean(1)


def mat_nosave(field):  # as mat but only final state (isolates save-buffer cost)
    Y0 = jnp.tile(y0[:, None], (1, n_v))
    def vf(t, Y, a):
        E = jnp.interp(t, ts, field)
        return C @ Y + E * (B @ Y) + (D @ Y) * deltas[None, :]
    return solve(vf, Y0, dfx.SaveAt(t1=True)).ys[0, 0].mean() * jnp.ones(n_t)


VARIANTS = {"base": base, "pfn": pfn, "mat": mat, "mat_nosave": mat_nosave}


def clock(f, *a):
    jax.block_until_ready(f(*a))  # warm
    c0, w0 = time.process_time(), time.perf_counter()
    jax.block_until_ready(f(*a))
    return time.perf_counter() - w0, time.process_time() - c0


name = sys.argv[1]
fn = VARIANTS[name]
fwd = jax.jit(fn)
grad = jax.jit(jax.grad(lambda f: jnp.sum(fn(f) ** 2)))
wf, cf = clock(fwd, field0)
wg, cg = clock(grad, field0)
print("RESULT " + json.dumps({
    "variant": name, "fwd_wall_s": wf, "fwd_cpu_s": cf, "grad_wall_s": wg, "grad_cpu_s": cg,
    "fwd_checksum": float(jnp.sum(fwd(field0))),
    "grad_checksum": float(jnp.sum(grad(field0))),
    "grad_norm": float(jnp.linalg.norm(grad(field0))),
}))
