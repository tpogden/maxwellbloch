"""Phase 3 headline numbers: forward-match errors, chain-gradient errors, reference timing."""
import json, sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "tests"); sys.path.insert(0, ".")
import numpy as np, jax, jax.numpy as jnp
import test_propagate as T
from mbjax.fd import central_fd, rel_err
from mbjax.lindblad import gaussian
from mbjax.propagate import build_propagation_model, propagate, velocity_classes

out = {}
# --- forward match vs clerq, velocity classes, 4 systems
for name, ampl, gg in T.CASES:
    spec = T._spec(name, ampl, detuning=0.2)
    kw = dict(t_min=-2.0, t_max=6.0, t_steps=160, z_steps=12, g=gg, vc=T.VC)
    ref = T._run_clerq(spec, **kw); got = np.asarray(T._run_jax(spec, ampl=ampl, **kw))
    out[f"fwd_{name}"] = {"rel_l2_real": rel_err(got.real, ref.real),
                          "ref_peak": float(np.abs(ref.real).max()),
                          "gain_vs_input": float(np.abs(ref.real[-1]).max() / np.abs(ref.real[0]).max())}
    print(name, out[f"fwd_{name}"], flush=True)
# --- chain gradients (two_level, 3 classes)
L, th = T._chain_loss("two_level", n_steps=T.N_FIXED)
g = np.asarray(jax.jit(jax.grad(L))(th)); fd = central_fd(L, th, list(range(len(th))), 1e-4)
Lp, _ = T._chain_loss("two_level")
gp = np.asarray(jax.jit(jax.grad(Lp))(th))
out["chain_grad"] = {"ad_vs_fd_fixed": rel_err(g, fd), "ad_pid_vs_ad_fixed": rel_err(gp, g)}
print(out["chain_grad"], flush=True)

# --- reference-size timing (two-level, n_z=200, n_t=1000, single class, AB, PID 1e-8)
spec = T._spec("two_level", 1.0)
model = build_propagation_model(spec)
ts = jnp.linspace(-3.0, 9.0, 1001); z_steps = 200; h = 1.0 / z_steps
N_z = jnp.ones(z_steps); deltas, w = velocity_classes(None)
pulse = gaussian(1.0, 0.0, 1.0)
def make(strategy):
    def loss(theta):
        om0 = 2 * np.pi * pulse(ts) * theta[0]
        om = propagate(model, om0, ts, N_z, 0.5 * theta[1], h, deltas, w, strategy=strategy,
                       rtol=1e-8, atol=1e-10)
        return jnp.sum(om[-1].real ** 2) * 0.012
    return loss
th0 = jnp.array([1.0, 1.0])
for strategy in ("s2", "s3"):
    L = make(strategy)
    f = jax.jit(L); vg = jax.jit(jax.value_and_grad(L))
    t0 = time.perf_counter(); jax.block_until_ready(f(th0)); c_f = time.perf_counter() - t0
    t0 = time.perf_counter(); jax.block_until_ready(vg(th0)); c_g = time.perf_counter() - t0
    ts_f, ts_g = [], []
    for _ in range(3):
        t0 = time.perf_counter(); jax.block_until_ready(f(th0)); ts_f.append(time.perf_counter() - t0)
        t0 = time.perf_counter(); jax.block_until_ready(vg(th0)); ts_g.append(time.perf_counter() - t0)
    out[f"timing_{strategy}"] = {"compile_fwd_s": c_f, "compile_grad_s": c_g,
                                 "fwd_s": min(ts_f), "grad_s": min(ts_g), "ratio": min(ts_g) / min(ts_f)}
    print(strategy, out[f"timing_{strategy}"], flush=True)
json.dump(out, open("bench/results/phase3_numbers.json", "w"), indent=1)
