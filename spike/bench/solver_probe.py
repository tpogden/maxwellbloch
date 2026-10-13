"""T4.4: accepted steps and time per slice for solver choices (one slice, forward only)."""
import json, sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import diffrax as dfx, jax, jax.numpy as jnp, numpy as np
from mbjax.lindblad import gaussian
from mbjax.propagate import build_propagation_model, velocity_classes, _rhs_batched
from mbjax.systems import make_system
from mbjax.chain import SYSTEM_PARAMS, T_MIN, T_MAX
jax.config.update("jax_enable_x64", True)

rows = []
for system in ("two_level", "j1_2_res"):
    ampl, g = SYSTEM_PARAMS[system]
    spec = make_system(system)
    spec["fields"][0].update(rabi_freq=1.0)
    model = build_propagation_model(spec)
    n_t, n_v = 1000, 10
    ts = jnp.linspace(T_MIN, T_MAX, n_t + 1)
    deltas, w = velocity_classes({"thermal_delta_min": -1.0, "thermal_delta_max": 1.0,
                                  "thermal_delta_steps": n_v - 1, "thermal_width": 0.7,
                                  "thermal_delta_inner_min": -1.0, "thermal_delta_inner_max": -1.0})
    pulse = gaussian(ampl, 0.0, 1.0)
    rho0 = jnp.broadcast_to(model["rho0"], (n_v,) + model["rho0"].shape)
    args = (model, lambda t: pulse(t), deltas, model["c"], model["cdc"], 1.0)
    for rtol in (1e-4, 1e-6, 1e-8):
        for name, solver in [("Tsit5", dfx.Tsit5()), ("Dopri5", dfx.Dopri5()), ("Dopri8", dfx.Dopri8()),
                             ("Bosh3", dfx.Bosh3())]:
            def run(r=rtol, s=solver):
                sol = dfx.diffeqsolve(dfx.ODETerm(_rhs_batched), s, T_MIN, T_MAX, 1e-3, rho0, args=args,
                                      saveat=dfx.SaveAt(t1=True),
                                      stepsize_controller=dfx.PIDController(rtol=r, atol=r * 1e-2),
                                      max_steps=200000)
                return sol.stats["num_steps"], sol.stats["num_rejected_steps"]
            try:
                f = jax.jit(run); jax.block_until_ready(f())
                t0 = time.perf_counter(); n, rej = [int(x) for x in jax.block_until_ready(f())]
                dt = time.perf_counter() - t0
                rows.append({"system": system, "rtol": rtol, "solver": name, "steps": n, "rejected": rej, "ms": dt * 1e3})
            except Exception as e:
                rows.append({"system": system, "rtol": rtol, "solver": name, "error": str(e)[:80]})
            print(rows[-1], flush=True)
json.dump(rows, open("bench/results/solver_probe.json", "w"), indent=1)
