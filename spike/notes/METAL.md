# Apple-silicon GPU (jax-metal): what was tried and what happened

Machine: Apple M1 Max, 64 GB. Scratch environment, separate from the spike's venv.

| Step | Result |
|-|-|
| `jax-metal` with jax / jaxlib 0.11.2 (spike's version) | Device detected; **every op fails** (`unknown attribute code: 22`), a StableHLO version mismatch with the plugin |
| `jax-metal` with jax / jaxlib 0.4.35 (the version it was built for) | float32 matmul, `while_loop`, `scan`, `grad` through `scan` **work**; **float64 and complex64 both fail** ("data types are not supported") |
| `diffrax` / `equinox` on jax 0.4.35 | **Cannot be imported**: every equinox that diffrax accepts needs a `jax.interpreters.batching.NotMapped` that jax 0.4.35 lacks; installing them upgrades jax and breaks the plugin again |
| Real-valued float32 RK4 Lindblad slice (`notes/metal_probe.py`: rho = A + iB, H real, sparse dissipator, 2400 steps, j1_2_res), forward | Same value as CPU f32 (0.26450 vs 0.26450 / 0.26697 vs 0.26697). 50 classes: **0.90 s Metal vs 0.55 s CPU**; 500 classes: **2.0 s Metal vs 4.7 s CPU** (2.3x faster) |
| Same, 2000 classes, forward | CPU 18.8 s; Metal **hung** (killed after > 10 min) |
| Same, **gradient** | **Fails on Metal**: `mps.multiply` with a null operand in the backward pass (reproduced standalone). A follow-up script trying variants (one-hot matmul instead of scatter, no `jax.checkpoint`) stalled for 10+ minutes at near-zero CPU before printing any result, so it is **not known** whether any variant works or which one stalled; not resolved |

Conclusion: the Apple GPU cannot run the differentiable backend in any form that matters
(f64 / complex are unsupported, `diffrax` is unusable with the only jax version the plugin
supports, and the backward pass fails or hangs). A float32 forward-only kernel can be faster
than the CPU for several hundred velocity classes, but float32 cannot meet the gradient-accuracy
gate, and the gradient is the point. A CUDA device (native float64) is the realistic GPU target.
