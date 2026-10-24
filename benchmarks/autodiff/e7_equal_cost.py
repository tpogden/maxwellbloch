"""E7: gradient against gradient-free search on the memory design problem, at equal cost.

All optimisers get the same budget, counted in units of one objective evaluation. A value-and-gradient
call is charged at the measured ratio (reverse-mode gradient time over forward time). Five optimisers:

  adam      optax Adam with a cosine learning-rate schedule, reverse-mode gradients
  lbfgsb    scipy L-BFGS-B on the reverse-mode gradient
  de        scipy differential evolution, one individual per parameter
  cmaes     CMA-ES (the `cma` package; run with `uv run --with cma ...`)
  lhs       Latin-hypercube random search

Each is started from, or sampled around, the same kind of random point (the start is the same for adam,
lbfgsb and cmaes at a given seed; de and lhs sample the box). The comparison follows the one in
Lee, Zenodo 10.5281/zenodo.23407812 (differential evolution and Latin hypercube against L-BFGS-B).

Run:  uv run --with cma python benchmarks/autodiff/e7_equal_cost.py --tag e7
"""

import argparse
import time

import numpy as np
from common import RESULTS, emit, time_once


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--evals-per-parameter", type=float, default=12.5)
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    p.add_argument(
        "--optimisers", nargs="+", default=["adam", "lbfgsb", "de", "cmaes", "lhs"]
    )
    p.add_argument(
        "--g",
        type=float,
        default=12.5,
        help="interaction strength; 12.5 is an optical depth of 25",
    )
    p.add_argument("--doppler", action="store_true")
    p.add_argument(
        "--box",
        type=float,
        default=4.0,
        help="half-width of the search box for de and lhs",
    )
    p.add_argument("--tag", default="e7")
    p.add_argument("--quick", action="store_true")
    a = p.parse_args()

    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    import optax
    from memory_problem import MemoryProblem
    from scipy.optimize import differential_evolution, minimize
    from scipy.stats import qmc

    out = RESULTS / f"{a.tag}.jsonl"
    prob = MemoryProblem(a.g, a.doppler)
    n = prob.n_par
    budget = 12.0 if a.quick else a.evals_per_parameter * n

    x0 = jnp.asarray(prob.random_start(0))
    time_once(lambda: prob.value(x0))
    time_once(lambda: prob.value_and_grad(x0))
    t_fwd = min(time_once(lambda: prob.value(x0)) for _ in range(3))
    t_vg = min(time_once(lambda: prob.value_and_grad(x0)) for _ in range(3))
    ratio = t_vg / t_fwd
    print(
        f"forward {t_fwd:.3f} s, value+grad {t_vg:.3f} s, charged {ratio:.2f} evaluations per gradient call",
        flush=True,
    )

    class Meter:
        """Counts cost in forward-evaluation units and records the best objective so far."""

        def __init__(self):
            self.cost, self.best, self.trace = 0.0, 0.0, [(0.0, 0.0)]

        def record(self, cost, value):
            self.cost += cost
            self.best = min(self.best, value)
            self.trace.append((self.cost, -self.best))

        def full(self):
            return self.cost >= budget

    def f(x, meter):
        v = float(prob.value(jnp.asarray(x)))
        meter.record(1.0, v)
        return v

    def fg(x, meter):
        v, g = prob.value_and_grad(jnp.asarray(x))
        meter.record(ratio, float(v))
        return float(v), np.asarray(g, dtype=float)

    class Stop(Exception):
        pass

    def run(name, seed):
        meter = Meter()
        start = prob.random_start(seed)
        rng_seed = 1000 + seed
        if name == "adam":
            steps = max(1, int(budget / ratio))
            opt = optax.adam(optax.cosine_decay_schedule(0.1, steps, alpha=0.05))
            coef, state = jnp.asarray(start), opt.init(jnp.asarray(start))
            for _ in range(steps):
                v, g = prob.value_and_grad(coef)
                meter.record(ratio, float(v))
                upd, state = opt.update(g, state)
                coef = optax.apply_updates(coef, upd)
            meter.record(1.0, float(prob.value(coef)))
        elif name == "lbfgsb":

            def fun(x):
                if meter.full():
                    raise Stop
                return fg(x, meter)

            try:
                minimize(
                    fun,
                    start,
                    jac=True,
                    method="L-BFGS-B",
                    options={"maxiter": 10_000, "maxfun": 10_000},
                )
            except Stop:
                pass
        elif name == "de":
            pop = (
                n  # popsize 1 means one individual per parameter in scipy's convention
            )
            generations = max(1, int(budget / pop) - 1)
            differential_evolution(
                lambda x: f(x, meter),
                [(-a.box, a.box)] * n,
                maxiter=generations,
                popsize=1,
                init="latinhypercube",
                seed=rng_seed,
                tol=0,
                polish=False,
                updating="immediate",
            )
        elif name == "cmaes":
            import cma

            es = cma.CMAEvolutionStrategy(
                start, 1.0, {"seed": rng_seed + 1, "verbose": -9}
            )
            while not es.stop() and meter.cost < budget:
                xs = es.ask()
                es.tell(xs, [f(x, meter) for x in xs])
        elif name == "lhs":
            pts = qmc.scale(
                qmc.LatinHypercube(d=n, seed=rng_seed).random(int(budget)),
                -a.box,
                a.box,
            )
            for x in pts:
                f(x, meter)
        else:
            raise ValueError(name)
        return meter

    for seed in a.seeds:
        for name in a.optimisers:
            t = time.perf_counter()
            meter = run(name, seed)
            step = max(1, len(meter.trace) // 60)
            emit(
                {
                    "exp": a.tag,
                    "optimiser": name,
                    "seed": seed,
                    "n_par": n,
                    "budget": budget,
                    "g": a.g,
                    "doppler": a.doppler,
                    "best_retrieved": meter.trace[-1][1],
                    "cost_used": meter.cost,
                    "charge_per_gradient": ratio,
                    "t_forward": t_fwd,
                    "t_value_and_grad": t_vg,
                    "wall_s": time.perf_counter() - t,
                    "trace": meter.trace[::step] + [meter.trace[-1]],
                },
                out,
            )
            print(
                f"seed {seed} {name:7s} best {max(r for _, r in meter.trace):.4f} cost {meter.cost:.0f} "
                f"wall {time.perf_counter() - t:.0f} s",
                flush=True,
            )


if __name__ == "__main__":
    main()
