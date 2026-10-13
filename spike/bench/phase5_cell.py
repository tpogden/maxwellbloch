"""Phase 5 cell: gradient cost of the reference chain for one (N, method).

methods: fwd (one forward), ad_s2 / ad_s3 (value_and_grad), fd_jax (central FD, 2N jitted
forwards), jacfwd (ForwardMode), fd_clerq (central FD through clerq PropagationSolver).
Reports best-of-``repeat`` wall-clock and the load average at start/end.
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import jax, jax.numpy as jnp, numpy as np
from mbjax.chain import make_reference_chain

jax.config.update("jax_enable_x64", True)
p = argparse.ArgumentParser()
p.add_argument("--n", type=int, required=True)
p.add_argument("--method", required=True)
p.add_argument("--repeat", type=int, default=3)
p.add_argument("--h", type=float, default=1e-5)
a = p.parse_args()
out = {"n": a.n, "method": a.method, "load_start": os.getloadavg()[0]}

def best(f, repeat):
    ts = []
    for _ in range(repeat):
        t0 = time.perf_counter(); r = f(); jax.block_until_ready(r); ts.append(time.perf_counter() - t0)
    return min(ts), ts

if a.method in ("fwd", "ad_s2", "ad_s3", "fd_jax", "jacfwd"):
    strategy = {"ad_s2": "s2", "ad_s3": "s3"}.get(a.method, "s2")
    loss, th, aux = make_reference_chain(a.n, strategy=strategy, forward_mode=(a.method == "jacfwd"))
    jax.block_until_ready(jax.jit(aux["output"])(th))  # builds the target outside timing
    _ = loss(th)
    if a.method == "fwd":
        f = jax.jit(loss); jax.block_until_ready(f(th))
        out["s"], out["s_all"] = best(lambda: f(th), a.repeat)
    elif a.method in ("ad_s2", "ad_s3"):
        f = jax.jit(jax.value_and_grad(loss))
        t0 = time.perf_counter(); jax.block_until_ready(f(th)); out["compile_run_s"] = time.perf_counter() - t0
        out["s"], out["s_all"] = best(lambda: f(th), a.repeat)
        out["grad_norm"] = float(jnp.linalg.norm(f(th)[1]))
    elif a.method == "jacfwd":
        f = jax.jit(jax.jacfwd(loss))
        t0 = time.perf_counter(); jax.block_until_ready(f(th)); out["compile_run_s"] = time.perf_counter() - t0
        out["s"], out["s_all"] = best(lambda: f(th), a.repeat)
        out["grad_norm"] = float(jnp.linalg.norm(f(th)))
    elif a.method == "fd_jax":
        f = jax.jit(loss); jax.block_until_ready(f(th))
        thn = np.asarray(th)
        def fd():
            g = np.empty(a.n)
            for i in range(a.n):
                e = np.zeros(a.n); e[i] = a.h
                g[i] = (float(f(jnp.asarray(thn + e))) - float(f(jnp.asarray(thn - e)))) / (2 * a.h)
            return jnp.asarray(g)
        out["s"], out["s_all"] = best(fd, max(1, min(a.repeat, 3 if a.n <= 50 else 1)))
        out["grad_norm"] = float(jnp.linalg.norm(fd()))
elif a.method == "fd_clerq":
    import copy
    from clerq import PropagationSolver
    from mbjax.systems import make_system
    from mbjax.chain import SYSTEM_PARAMS, T_MIN, T_MAX
    from mbjax.lindblad import gaussian
    loss, th, aux = make_reference_chain(a.n)
    ampl, g = SYSTEM_PARAMS["two_level"]
    spec = make_system("two_level")
    spec["fields"][0].update(rabi_freq=1.0, rabi_freq_t_func="gaussian",
                             rabi_freq_t_args={"ampl": ampl, "centre": 0.0, "fwhm": 1.0})
    basis = np.asarray(aux["basis"]); ts = np.asarray(aux["ts"])
    centres_np = np.asarray(aux["centres"]); sigma = aux["sigma"]
    target = np.asarray(jax.jit(aux["output"])(aux["target_theta"]))
    dt = (T_MAX - T_MIN) / 1000
    def clerq_loss(theta):
        ps = PropagationSolver(system=copy.deepcopy(spec), t_min=T_MIN, t_max=T_MAX, t_steps=1000,
                               z_min=0.0, z_max=1.0, z_steps=200, z_steps_inner=1,
                               interaction_strengths=[g],
                               opts={"atol": 1e-8, "rtol": 1e-6, "method": "adams", "nsteps": 100000})
        th_np = np.asarray(theta)
        # give clerq the pulse as its analytic t_func (cyclic units; rabi_freq = 1)
        f0 = ps.system.fields[0]
        f0.rabi_freq = 1.0
        f0.rabi_freq_t_func = lambda t, args=None: (
            np.exp(-((np.asarray(t)[..., None] - centres_np) ** 2) / (2 * sigma**2)) @ th_np
            if np.ndim(t) else float(np.exp(-((t - centres_np) ** 2) / (2 * sigma**2)) @ th_np))
        ps.system._build_H_Omega()
        ps.solve(step="ab", recalc=True, progress=False)
        return float(np.sum((ps.Omegas_zt[0, -1].real / (2 * np.pi) - target) ** 2) * dt), ps
    # one forward (the unit cost), and the equivalence check vs the JAX loss
    ts_f = []
    for _ in range(2):
        t0 = time.perf_counter(); lv, _ps = clerq_loss(np.asarray(th)); ts_f.append(time.perf_counter() - t0)
    out["fwd_s"] = min(ts_f)
    out["loss_clerq"] = lv; out["loss_jax"] = float(loss(th))
    out["s"] = 2 * a.n * min(ts_f)          # extrapolated: 2N forwards
    out["extrapolated"] = True
    if a.n <= 100:                           # up to N=100: run the real FD
        thn = np.asarray(th); t0 = time.perf_counter()
        gg = np.empty(a.n)
        for i in range(a.n):
            e = np.zeros(a.n); e[i] = a.h
            gg[i] = (clerq_loss(thn + e)[0] - clerq_loss(thn - e)[0]) / (2 * a.h)
        out["s_measured"] = time.perf_counter() - t0
        out["extrapolated"] = False
        out["s"] = out["s_measured"]
        out["grad_norm"] = float(np.linalg.norm(gg))
out["load_end"] = os.getloadavg()[0]
print("RESULT " + json.dumps(out))
