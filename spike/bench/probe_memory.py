"""T0.2: does jax.checkpoint on the scan body change memory vs plain scan?

Toy chain: for each z-slice, solve y' = (-g I + E(t) B + C) y for a dim^2-vector y,
driven by a linearly interpolated field E(z, t); polarisation P(t) = y[0]; update
E <- E + dz * kappa * P. Loss = sum(E_out^2).

Strategies (inner is always Tsit5 + ConstantStepSize with n_t steps):
  direct    : DirectAdjoint (stores every step), plain scan
  s2        : RecursiveCheckpointAdjoint, plain scan
  s3        : RecursiveCheckpointAdjoint, jax.checkpoint on scan body
  s5        : s3 + chunked (nested) scan over z, chunk ~ sqrt(n_z)

Reports compile-time temp bytes from XLA's memory_analysis and, with --rss,
peak RSS of a fresh subprocess running the compiled value_and_grad.
"""

import argparse
import json
import math
import resource
import subprocess
import sys
import time

import diffrax as dfx
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)


def make_problem(n_z, n_t, dim, strategy, n_v=1, seed=0):
    n = dim * dim
    kB, kC = jax.random.split(jax.random.PRNGKey(seed))
    B = 0.3 * jax.random.normal(kB, (n, n))
    C = 0.3 * (jax.random.normal(kC, (n, n)))
    C = C - C.T - 0.5 * jnp.eye(n)  # antisym + decay: stable
    ts = jnp.linspace(0.0, 1.0, n_t)
    y0 = jnp.zeros(n).at[0].set(1.0)
    dz, kappa = 1.0 / n_z, 0.5

    adjoint = (
        dfx.DirectAdjoint()
        if strategy == "direct"
        else dfx.RecursiveCheckpointAdjoint()
    )

    deltas = jnp.linspace(-1.0, 1.0, n_v)  # per-velocity-class detuning
    D = 0.3 * jax.random.normal(jax.random.PRNGKey(seed + 1), (n, n))
    D = D - D.T

    def solve_class(field, delta):
        def vf(t, y, args):
            E = jnp.interp(t, ts, field)
            return (C + E * B + delta * D) @ y

        sol = dfx.diffeqsolve(
            dfx.ODETerm(vf),
            dfx.Tsit5(),
            0.0,
            1.0,
            1.0 / (n_t - 1),
            y0,
            saveat=dfx.SaveAt(ts=ts),
            stepsize_controller=dfx.ConstantStepSize(),
            adjoint=adjoint,
            max_steps=n_t + 8,
        )
        return sol.ys[:, 0]  # polarisation(t)

    def solve_slice(field):
        # vmap over velocity classes inside the slice body; carry stays n_t
        return jax.vmap(solve_class, in_axes=(None, 0))(field, deltas).mean(0)

    def body(field, _):
        P = solve_slice(field)
        return field + dz * kappa * P, None

    if strategy in ("s3", "s5"):
        body = jax.checkpoint(body)

    def propagate(field0):
        if strategy == "s5":
            chunk = max(1, int(round(math.sqrt(n_z))))
            n_chunks, rem = divmod(n_z, chunk)

            @jax.checkpoint
            def chunk_body(field, _):
                field, _ = jax.lax.scan(body, field, None, length=chunk)
                return field, None

            field, _ = jax.lax.scan(chunk_body, field0, None, length=n_chunks)
            if rem:
                field, _ = jax.lax.scan(body, field, None, length=rem)
            return field
        field, _ = jax.lax.scan(body, field0, None, length=n_z)
        return field

    def loss(field0):
        return jnp.sum(propagate(field0) ** 2)

    field0 = jnp.exp(-((ts - 0.5) ** 2) / 0.02)
    return jax.value_and_grad(loss), field0


def run_cell(args):
    vg, field0 = make_problem(args.n_z, args.n_t, args.dim, args.strategy, args.n_v)
    vg = jax.jit(vg)
    t0 = time.perf_counter()
    compiled = vg.lower(field0).compile()
    t_compile = time.perf_counter() - t0
    ma = compiled.memory_analysis()
    out = {
        "strategy": args.strategy,
        "n_z": args.n_z,
        "n_t": args.n_t,
        "dim": args.dim,
        "n_v": args.n_v,
        "temp_MiB": ma.temp_size_in_bytes / 2**20,
        "arg_MiB": ma.argument_size_in_bytes / 2**20,
        "compile_s": t_compile,
    }
    if args.rss:
        v, g = compiled(field0)
        jax.block_until_ready(g)
        t0 = time.perf_counter()
        v, g = compiled(field0)
        jax.block_until_ready(g)
        out["grad_s"] = time.perf_counter() - t0
        out["loss"] = float(v)
        out["grad_norm"] = float(jnp.linalg.norm(g))
        # ru_maxrss is bytes on macOS, KiB on Linux
        scale = 1 if sys.platform == "darwin" else 1024
        out["peak_rss_MiB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 2**20
    print("RESULT " + json.dumps(out))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n_z", type=int, default=100)
    p.add_argument("--n_t", type=int, default=500)
    p.add_argument("--dim", type=int, default=4)
    p.add_argument("--n_v", type=int, default=1)
    p.add_argument("--strategy", default="s2", choices=["direct", "s2", "s3", "s5"])
    p.add_argument("--rss", action="store_true")
    p.add_argument("--sweep", action="store_true", help="run a grid in subprocesses")
    p.add_argument("--out", default="bench/results/probe_memory.json")
    args = p.parse_args()

    if not args.sweep:
        return run_cell(args)

    rows = []
    grid = [(nz, 500, 4) for nz in (50, 100, 200, 400)] + [(100, 500, d) for d in (4, 9)]
    for n_z, n_t, dim in grid:
        for strat in ("direct", "s2", "s3", "s5"):
            cmd = [sys.executable, __file__, "--n_z", str(n_z), "--n_t", str(n_t),
                   "--dim", str(dim), "--strategy", strat, "--rss"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
            if r.returncode or not line:
                rows.append({"strategy": strat, "n_z": n_z, "n_t": n_t, "dim": dim,
                             "error": r.stderr[-300:]})
            else:
                rows.append(json.loads(line[0][7:]))
            print(rows[-1], flush=True)
    with open(args.out, "w") as f:
        json.dump(rows, f, indent=1)


if __name__ == "__main__":
    main()
