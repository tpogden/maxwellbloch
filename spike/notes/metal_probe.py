"""Fixed-step RK4 Lindblad slice, real-valued float32 (rho = A + iB, H real): CPU vs Apple GPU."""
import sys, time, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, jax, jax.numpy as jnp
ops = np.load(sys.argv[3])
f32 = jnp.float32
n = ops["X"].shape[0]
X, N = jnp.asarray(ops["X"], f32), jnp.asarray(ops["N"], f32)
a, b, w, r = jnp.asarray(ops["a"]), jnp.asarray(ops["b"]), jnp.asarray(ops["w"], f32), jnp.asarray(ops["r"], f32)
rows, cols, cw = jnp.asarray(ops["rows"]), jnp.asarray(ops["cols"]), jnp.asarray(ops["cw"], f32)
PI = np.pi
T0, T1, N_OUT, INNER = -3.0, 9.0, 240, 10            # 2400 RK4 steps, dt = 0.005
DT = (T1 - T0) / (N_OUT * INNER)

def make(n_v):
    deltas = jnp.linspace(-1.0, 1.0, n_v).astype(f32)
    wts = jnp.ones(n_v, f32) / n_v
    eye = jnp.eye(n, dtype=f32)
    damp = -0.5 * (r[:, None] + r[None, :])
    def rhs(t, y, ampl):
        A, B = y[0], y[1]
        om = ampl * jnp.exp(-4.0 * jnp.log(2.0) * t * t)
        H = PI * om * X[None] - 2 * PI * deltas[:, None, None] * N[None]
        comm = lambda M: H @ M - M @ H
        gain = jnp.zeros((n_v, n), f32).at[:, a].add(w * A[:, b, b])
        return jnp.stack([comm(B) + gain[:, :, None] * eye + damp * A, -comm(A) + damp * B])
    def step(y, t, ampl):
        k1 = rhs(t, y, ampl); k2 = rhs(t + DT / 2, y + DT / 2 * k1, ampl)
        k3 = rhs(t + DT / 2, y + DT / 2 * k2, ampl); k4 = rhs(t + DT, y + DT * k3, ampl)
        return y + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    A0 = jnp.zeros((n_v, n, n), f32).at[:, jnp.arange(3), jnp.arange(3)].set(1.0 / 3)
    y0 = jnp.stack([A0, jnp.zeros_like(A0)])
    def loss(ampl):
        @jax.checkpoint
        def chunk(y, i):
            t0 = T0 + i * INNER * DT
            y = jax.lax.fori_loop(0, INNER, lambda j, yy: step(yy, t0 + j * DT, ampl), y)
            coh = jnp.einsum("v,vk->", wts, y[1][:, rows, cols] * cw[None])   # Im part: absorption
            return y, coh
        _, cs = jax.lax.scan(chunk, y0, jnp.arange(N_OUT, dtype=f32))
        return jnp.sum(cs) * (T1 - T0) / N_OUT
    return loss

def timeit(f, x, rep=3):
    jax.block_until_ready(f(x)); best = 1e9
    for _ in range(rep):
        t0 = time.perf_counter(); jax.block_until_ready(f(x)); best = min(best, time.perf_counter() - t0)
    return best

for backend in sys.argv[1].split(","):
    dev = jax.devices("cpu")[0] if backend == "cpu" else [d for d in jax.devices() if d.platform.lower() == "metal"][0]
    for n_v in [int(v) for v in sys.argv[2].split(",")]:
        row = {"backend": backend, "n_v": n_v}
        try:
            with jax.default_device(dev):
                loss = make(n_v); x = jnp.asarray(1.5, f32)
                f = jax.jit(loss); row["value"] = float(f(x)); row["fwd_s"] = round(timeit(f, x), 3)
                g = jax.jit(jax.value_and_grad(loss)); v, gr = g(x); row["grad"] = float(gr); row["grad_s"] = round(timeit(g, x), 3)
        except Exception as e:
            row["error"] = f"{type(e).__name__}: {str(e).splitlines()[0][:150]}"
        print(json.dumps(row), flush=True)
