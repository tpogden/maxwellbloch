"""T5.2: inverse design on the reference medium: recover an input pulse that produces a target
output field, with L-BFGS-B driven by (a) AD gradients, (b) JAX finite differences; (c) clerq-FD
cost is projected from the measured clerq forward time and the number of gradient evaluations.

Planted target: the output of a known smooth input theta_true, so zero loss is attainable.
"""
import json, sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import jax, jax.numpy as jnp, numpy as np
from scipy.optimize import minimize
from mbjax.chain import make_reference_chain

jax.config.update("jax_enable_x64", True)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
T_CLERQ_FWD = float(sys.argv[2]) if len(sys.argv) > 2 else 2.25   # measured clerq forward (s)
H_FD = 1e-5

_, th0, aux0 = make_reference_chain(N)
k = np.arange(N)
theta_true = np.asarray(th0) * (1.0 + 0.7 * np.cos(2 * np.pi * 1.5 * k / N) * np.exp(-((k - N / 2) / (N / 3)) ** 2))
theta_true = jnp.asarray(theta_true)
loss, th0, aux = make_reference_chain(N, target_theta=theta_true)
loss_j = jax.jit(loss)
vg = jax.jit(jax.value_and_grad(loss))
jax.block_until_ready(vg(th0)); jax.block_until_ready(loss_j(th0))     # compile outside timing
L0 = float(loss_j(th0))
print(f"N={N}  initial loss {L0:.3e}", flush=True)

def run(kind, maxiter=80, gtol=1e-12):
    hist = []
    t_start = time.perf_counter()
    n_grad = [0]

    def fun(x):
        x = jnp.asarray(x)
        if kind == "ad":
            v, g = vg(x)
            g = np.asarray(g)
        else:
            v = float(loss_j(x))
            xn = np.asarray(x); g = np.empty(N)
            for i in range(N):
                e = np.zeros(N); e[i] = H_FD
                g[i] = (float(loss_j(jnp.asarray(xn + e))) - float(loss_j(jnp.asarray(xn - e)))) / (2 * H_FD)
        n_grad[0] += 1
        hist.append((time.perf_counter() - t_start, float(v)))
        return float(v), g

    res = minimize(fun, np.asarray(th0), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15, "gtol": gtol, "maxfun": 200})
    return {"kind": kind, "loss": float(res.fun), "nit": int(res.nit), "n_grad": n_grad[0],
            "wall_s": hist[-1][0], "hist": hist, "x": res.x.tolist(), "message": str(res.message)}

out = {"N": N, "L0": L0, "runs": {}}
for kind in ("ad", "fd_jax"):
    r = run(kind)
    out["runs"][kind] = r
    print(f"{kind:7} nit {r['nit']:3d} grad evals {r['n_grad']:3d}  loss {r['loss']:.3e} "
          f"({L0 / max(r['loss'], 1e-300):.1e}x lower)  wall {r['wall_s']:.1f}s", flush=True)
ad = out["runs"]["ad"]
out["projected_clerq_fd_s"] = ad["n_grad"] * 2 * N * T_CLERQ_FWD
print(f"clerq-FD projected: {ad['n_grad']} gradients x {2*N} forwards x {T_CLERQ_FWD:.2f}s = "
      f"{out['projected_clerq_fd_s']/60:.1f} min", flush=True)
np.savez("bench/results/optimise.npz", theta_true=np.asarray(theta_true), theta0=np.asarray(th0),
         theta_ad=np.asarray(ad["x"]), theta_fd=np.asarray(out["runs"]["fd_jax"]["x"]),
         out_true=np.asarray(jax.jit(aux["output"])(theta_true)),
         out_init=np.asarray(jax.jit(aux["output"])(th0)),
         out_ad=np.asarray(jax.jit(aux["output"])(jnp.asarray(ad["x"]))),
         ts=np.asarray(aux["ts"]))
json.dump(out, open("bench/results/optimise.json", "w"), indent=1)
