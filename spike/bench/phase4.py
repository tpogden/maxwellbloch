"""Phase 4 cell runner: one (system, n_z, n_t, n_v, strategy) -> JSON line.

Each cell is meant to run in its own subprocess (clean peak RSS). Reports compile-time temp
memory (XLA memory_analysis), peak RSS, forward / value_and_grad wall-clock (best of
``--repeat`` after warm-up), and accepted ODE steps for slice 0 of the forward pass.
"""
import argparse, json, resource, sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import jax, jax.numpy as jnp
from mbjax.chain import make_chain

jax.config.update("jax_enable_x64", True)
p = argparse.ArgumentParser()
p.add_argument("--system", default="j1_2_res")
p.add_argument("--n_z", type=int, default=200)
p.add_argument("--n_t", type=int, default=1000)
p.add_argument("--n_v", type=int, default=10)
p.add_argument("--strategy", default="s3", choices=["s1", "s2", "s3", "s4", "s5", "fwd", "direct"])
p.add_argument("--checkpoints", type=int, default=None)
p.add_argument("--chunk", type=int, default=None)
p.add_argument("--max_steps", type=int, default=16384)
p.add_argument("--rtol", type=float, default=1e-6)
p.add_argument("--repeat", type=int, default=2)
p.add_argument("--compile_only", action="store_true")
p.add_argument("--no_warmup", action="store_true", help="time the first execution (long runs)")
a = p.parse_args()

strategy = a.strategy
ckpt = a.checkpoints if strategy == "s4" else None
if strategy == "s4":
    strategy_impl = "s3"
else:
    strategy_impl = strategy
loss, th = make_chain(a.system, a.n_z, a.n_t, a.n_v, strategy=strategy_impl, rtol=a.rtol,
                      atol=a.rtol * 1e-2, max_steps=a.max_steps, checkpoints=ckpt,
                      chunk=a.chunk, forward_mode=(strategy == "fwd"))
out = {k: getattr(a, k) for k in ("system", "n_z", "n_t", "n_v", "strategy", "max_steps", "rtol",
                                  "checkpoints", "chunk")}

def fn_for(strategy):
    if strategy == "s1":
        return jax.jit(loss)
    if strategy == "fwd":
        return jax.jit(jax.jacfwd(loss))
    return jax.jit(jax.value_and_grad(loss))

f = fn_for(strategy)
t0 = time.perf_counter()
c = f.lower(th).compile()
out["compile_s"] = time.perf_counter() - t0
ma = c.memory_analysis()
out["temp_MiB"] = ma.temp_size_in_bytes / 2**20
if not a.compile_only:
    if not a.no_warmup:
        jax.block_until_ready(c(th))
    ts = []
    for _ in range(a.repeat):
        t0 = time.perf_counter(); jax.block_until_ready(c(th)); ts.append(time.perf_counter() - t0)
    out["run_s"] = min(ts)
    out["run_s_all"] = ts
    scale = 1 if sys.platform == "darwin" else 1024
    out["peak_rss_MiB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 2**20
    if strategy != "s1":
        g = c(th)
        out["grad_norm"] = float(jnp.linalg.norm(g[1] if isinstance(g, tuple) else g))
print("RESULT " + json.dumps(out))
