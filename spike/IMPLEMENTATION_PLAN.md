# Implementation plan: a differentiable propagation backend for clerq

Follows from the spike ([RESULTS.md](RESULTS.md), decision: build). This plan is a proposal for
review; items marked **DECISION** need an answer before work starts. Effort figures are
rough estimates for one person, not measurements.

## 1. Goal and scope

**Goal.** Let users compute exact gradients of a scalar loss of the propagated field with
respect to problem parameters (input pulse, interaction strength, decay rates, detuning), by
reverse-mode autodiff through the same Maxwell-Bloch scheme `PropagationSolver` uses, so that
gradient-based design of pulses and media is tractable.

**First cut (v1):**

- Co-propagating propagation with the existing finite-difference scheme (Adams-Bashforth
  default, Euler available) and velocity classes, for systems and fields as described by clerq's
  `System` / `Field` objects and JSON problem files.
- float64 on CPU (CUDA float64 is the intended GPU path; not yet validated).
- Reverse mode (default) and forward mode.
- Results are validated against `PropagationSolver` as a live oracle.

**Non-goals for v1:** counter-propagating fields (needs a boundary-value iteration, not a scan),
non-co-propagating geometries, replacing `PropagationSolver`, the Apple GPU (see
`notes/METAL.md`: no float64/complex support and the plugin cannot run our stack), float32 as a
supported precision.

## 2. Decisions needed

| # | Decision | Options | Recommendation |
|-|-|-|-|
| D1 | API shape | (a) new `AutodiffSolver` class repeating the constructor args; (b) `PropagationSolver.to_autodiff()` returning a pure-function object built from the same problem definition | **Decided: (b).** One problem definition, JSON files work unchanged, no argument-list drift |
| D2 | Where it lives | `clerq.autodiff` subpackage (D6); `jax`, `diffrax`, `equinox` as an optional extra `clerq[jax]` | **Decided: yes.** Never imported by default, so the core install stays QuTiP + NumPy + SciPy |
| D3 | Number of fields in v1 | one / several co-propagating | **Decided: several co-propagating.** The spike's operators already carry a field axis; the work is the per-field source and tests (folded into M3) |
| D4 | float64 policy | library enables `jax_enable_x64` globally on first use / requires the user to / raises if disabled | **Decided: raise** with a clear message if x64 is off; never mutate global JAX config silently |
| D5 | Python floor | clerq declares >= 3.10 | **Verify** that the pinned JAX/diffrax support it; if not, the extra declares a higher floor |
| D6 | Name | `clerq.autodiff` (chosen) | **Decided: `clerq.autodiff`** |

## 3. Design

### What is ported from the spike, and what is left behind

| Spike module | Becomes | Change needed |
|-|-|-|
| `mbjax/systems.py::operators` | `clerq/autodiff/_operators.py` | Build X_k, N_k and collapse operators **from a `System` object** (not a spec dict); keep the single-element detection for the sparse dissipator, dense fallback |
| `mbjax/lindblad.py` | `_lindblad.py` | Drop spec/clerq-import coupling; keep `check_fixed_step` |
| `mbjax/drive.py` | `_drive.py` | Keep the not-a-knot spline; replace dense derivative operator by a banded solve before n_t grows |
| `mbjax/propagate.py` | `_propagate.py` | Generalise from field 0 to a field axis; keep AB2/Euler, velocity batching, S2/S3/S4 |
| `mbjax/chain.py`, `bench/` | `tests/bench_autodiff.py` | Benchmarks only |
| `tests/` | `tests/test_autodiff.py` | Port, with the live-oracle design below |

Left behind: the system specs (clerq has its own), the `legacy` Doppler / source-weight flags
(master has the corrected conventions; carry the legacy ones only if someone needs old results),
S5 chunked scan unless a use case appears (S3/S4 sufficed everywhere we measured).

### Data flow

```
JSON / dict problem ──> PropagationSolver ──(.to_autodiff())──> AutodiffPropagation
                                                                      │
   params pytree {pulse, g_scale, gamma_scale, detuning, ...} ───────>│ propagate(params) -> Omega(z,t)
                                                                      │ value_and_grad(loss, params)
```

- `params` is a plain pytree of JAX arrays; the object closes over everything static (operators,
  grids, velocity classes, number density profile).
- `propagate` is pure and `jit`-able, so users can compose it with `optax`, `jax.vmap` (batch of
  problems) or their own losses.
- Defaults from the spike: Tsit5 + PID control (rtol 1e-6 forward, tighter for gradient checks),
  S3 checkpointed scan body, `max_steps` derived from the pulse and exposed, return only the last
  slice unless `return_all=True`.

### Conventions the port must reproduce (all verified in the spike against clerq master)

- Only Re(Omega) drives H; Omega itself is complex.
- Slice 0 uses the analytic pulse, later slices a cubic **not-a-knot** spline of the sampled field.
- Source term weights each coupled coherence by its signed factor; N taken at the next z in both
  Adams-Bashforth terms; one Euler step bootstraps AB.
- Velocity-class detunings are cyclic `thermal_delta / (2 pi)`.
- `states_zt` is **not** reproduced (off by one in clerq); only `Omegas_zt` is a contract.

Because these have already changed twice, the test suite uses `PropagationSolver` as a live
oracle (section 5): a future clerq convention change fails the tests instead of silently
diverging.

## 4. Milestones

| M | Deliverable | Acceptance | Est. |
|-|-|-|-|
| M0 | Scaffolding: `clerq[jax]` extra, `clerq.autodiff` skeleton, version pins and Python floor (D5), CI job for the extra, x64 policy (D4) | `pip install clerq[jax]` works on CI; core install and tests unaffected | 1 d |
| M1 | Operators from `System` | X_k, N_k, collapse operators equal `System.H_Omega_list`, `H_Delta`, `c_ops` to 1e-12 for clerq's fixtures | 1 d |
| M2 | Slice solve + array drive | Matches `LindbladSolver` to 1e-8 with analytic and interpolated drive; gradients w.r.t. field samples vs fixed-step FD <= 1e-5 | 2 d |
| M3 | Propagation scan (AB2/Euler, velocity classes, **several co-propagating fields**, strategies) | Matches `PropagationSolver` on the fixture set in `tests/json/` to 1e-8 relative (same scheme, same z_steps), including the multi-field fixtures | 4-5 d |
| M4 | Differentiation API, parameter pytree, forward mode | Gradients vs fixed-step FD <= 1e-5; reverse vs forward <= 1e-8; strategies agree to 1e-9 | 2-3 d |
| M5 | Hardening: guards (commensurate fixed step, non-vacuity), `max_steps` heuristic, clear errors for unsupported features (counter-propagating, ...), memory/time docs | Guard tests; unsupported inputs raise `NotImplementedError` with guidance | 2 d |
| M6 | Docs and example: inverse-design notebook, API page, CHANGELOG + release notes, benchmark file | `make docs_html` clean; example runs in CI time budget | 2-3 d |
| M7 (optional) | CUDA float64 validation and timing | One sweep on a GPU; documented speedups | 1 d + hardware |

Rough total for M0-M6: 2.5-3.5 weeks.

### Status and M5 checklist

M0-M4 are **merged to master** in PR tpogden/clerq#327 (merge commit 3c46011, 2026-10-15),
including the CHANGELOG and release-notes entries that were planned for M6. The clerq
Adams-Bashforth slice-offset defect found while porting was fixed separately in #328 and
the tests now compare against `PropagationSolver` directly for every `z_steps_inner`. A code
review of the PR (10 findings) produced the hardening list below; "done" items were fixed
before the merge, the rest are M5 unless noted. M5 work continues on `feature/autodiff-m5`.

| # | Item | Status | Notes |
|-|-|-|-|
| R1 | A `square` pulse's `on`/`off` have silently-zero gradients | **done** | `value_and_grad` holds them fixed by default, raises if they are requested, takes `fixed=...`; ramp pulses stay differentiable |
| R2 | Integer literals in a problem definition break `value_and_grad(loss)` | **done** | `default_params()` returns floats; integer user parameters are converted |
| R3 | `value_and_grad` compile cache is unbounded and keyed by the loss closure | M5 | Bound it (small LRU) or key on something stable; document reuse of the loss object |
| R4 | Dense n_t x n_t spline operator (about 200 MB at n_t = 5000, quadratic) | M5 | Banded / tridiagonal-plus-corrections solve; already listed in the risk table |
| R5 | Pulse ports re-implement `t_funcs`; `ramp_on`, `ramp_off`, `ramp_offon` and `sech(fwhm)` have no direct equivalence test | M5 | Add per-function tests against `clerq.t_funcs` on a time grid (or make `t_funcs` backend-agnostic) |
| R6 | `strategy="s4"` accepted without `checkpoints`; default `max_steps=100_000` unrelated to memory | M5 | Validate `checkpoints` for s4; derive a `max_steps` default from the pulse and document the s2 memory coupling |
| R7 | `omega0` is angular while every other parameter is cyclic | **done** | Decided: cyclic, converted internally (x 2 pi); the returned field keeps the units of `Omegas_zt` (angular), stated in the docstring; guard test for the 2 pi error |
| R8 | Problem with no fields crashes with an obscure `np.concatenate` error | M5 | Raise a clear message (or support the empty case) |
| R9 | Pulse-argument exclusivity (ampl/n_pi, width/fwhm) encoded in two modules | M5, low | Declare the alternative groups once next to the pulse definitions |
| R10 | Broken `:class:`cls`` docstring reference | M6 | Resolved when `clerq.autodiff` is added to the API docs |
| R11 | True derivative w.r.t. a square pulse's edge times (about -3590 / +3058 on the example) | optional, after M5 | A feature rather than hardening: integrate slice 0 in segments with differentiable endpoints (diffrax differentiates integration limits); the sampled input field is piecewise constant in the edge times, so only slice 0 contributes |

Findings that are guards, not behaviour, are covered by tests in the PR
(`tests/test_autodiff_api_regressions.py` for R1 and R2).

## 5. Test plan

- **Live oracle.** Each supported fixture in `tests/json/` is solved with `PropagationSolver`
  at test time and compared (no stored golden arrays, so no stale references).
- **Non-vacuity guard.** Every oracle comparison asserts the medium changes the field by more
  than 1% in the oracle itself (the spike nearly shipped a "matches to 3e-16" test where both
  sides returned zero).
- **Gradient checks.** Central finite differences on a **fixed-step** solve with a step count
  that is not commensurate with the output grid; adaptive-solver gradients are compared with the
  fixed-step AD, never with adaptive FD (it is noise-limited).
- **Strategy and mode agreement.** S2 / S3 / S4 and reverse / forward agree to 1e-9 / 1e-8.
- **Guards.** The commensurate-step check raises; `x64` disabled raises.
- **Markers and runtime.** `@pytest.mark.jax`, skipped when the extra is absent; problems sized
  so the whole file stays within a few minutes (the spike suite is ~7 min at its sizes and can
  be pruned).
- **Benchmarks** in `tests/bench_autodiff.py` (the repo already has `bench_*.py`), not
  part of the default run.

## 6. Risks and mitigations

| Risk | Evidence | Mitigation |
|-|-|-|
| JAX / diffrax API churn and version coupling | `jax-metal` forced jax 0.4.35 where `diffrax` could not even be imported; the spike pinned jax 0.11.2 / diffrax 0.7.2 | Pin a tested range in the extra; run CI against the lowest and highest; keep the surface that touches diffrax small |
| A diffrax gradient bug on commensurate fixed steps | `notes/diffrax_commensurate_repro.py` | Guard in code; report upstream; test that the guard fires |
| clerq conventions change again | two defects found and fixed during the spike | Live-oracle tests fail loudly; one place documents the contract |
| First-call compile time | ~1-5 s per new shape | Document; `jit` once; keep shapes static |
| Dense spline operator at large n_t | 200 MB at n_t = 5000 | Banded solve in M2/M5 |
| Memory under adaptive stepping with plain scan | memory scales with `max_steps` | Default to the checkpointed body; derive and expose `max_steps` |
| Scope creep to counter-propagation / dense dissipators / stiff regimes | all explicitly unmeasured in the spike | Non-goals for v1; raise clear errors |
| GPU expectations | No CUDA numbers; Apple GPU not viable | State plainly in docs; M8 measures it |

## 7. Open questions

1. D1-D6 above.
2. Which clerq fixtures should be the oracle set (all of `tests/json/`, or a curated subset that
   runs quickly)?
3. Should the backend accept the **same** `rabi_freq_t_func` names (`gaussian`, `sech`,
   `square`) as JAX-traceable functions so that pulse parameters are differentiable, or only
   sampled arrays in v1? (The spike used sampled arrays plus an analytic slice-0 pulse; making
   `t_funcs` traceable is more work but gives users parameters like width and amplitude.)
4. How do we want to expose the number-density profile (`num_density_z_func`) as a
   differentiable parameter, if at all, in v1?
5. Who validates the CUDA path, and on what hardware?
