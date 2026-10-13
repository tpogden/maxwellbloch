# Spike task plan: differentiable Maxwell-Bloch backend

Executes [DESIGN.md](DESIGN.md) with the amendments below. Go/no-go on a JAX/diffrax
backend for clerq. The package itself is not modified by this spike.

## Decisions and amendments to DESIGN.md

| # | Amendment | Why |
|---|-----------|-----|
| A1 | Run in float64/complex128 on **CPU** (M1 Max, 64 GB RAM). Metal has no f64. | Gates 1-3 need f64. |
| A2 | Gate 1 reworded: S4 fits the stretch point in <= 64 GB host RAM, projected <= 32 GB device memory. CUDA timing run deferred until gates pass. | No discrete GPU locally. |
| A3 | Gate 2 replaced: peak memory at n_z=2000 <= 4x peak at n_z=200, measured at stretch-like dim/n_v, fit restricted to the regime above the fixed memory floor. | Carry is the field (n_t), so alpha ~ 1 with a tiny constant is fine; alpha alone can fail for the wrong reason. |
| A4 | Add strategies **S5** (chunked/nested scan, sqrt(n_z)) and **FWD** (`jax.jacfwd`). | Plain `jax.checkpoint` on the body is O(n_z) in carries; S2 and S3 may be near-identical. |
| A5 | Gate 3: AD vs fixed-step, tight-tolerance central FD, and AD vs `jacfwd`. Rel. error <= 1e-5. | Adaptive steps make FD noisy. |
| A6 | Gate 4 also reports AD vs `jacfwd` at N=10-100; add wall-clock gate: backward/forward <= 10x at reference. | Gate 4 as written is ~arithmetic (2N/c). |
| A7 | Doppler classes: `vmap` **inside** the slice body, sum P over v, carry stays n_t. Compare adaptive vs fixed-step under vmap. | Keeps n_v multiplier on recomputed, not stored, memory. |
| A8 | Reference problem fixed (below). n_params = 50 pulse-basis coefficients; full n_t-sample field gradient is a separate case. | DESIGN.md reference was internally inconsistent. |
| A9 | Fields as real/imag pairs. | Avoids JAX complex-gradient conjugate convention. |

### Reference problem (single source of truth)

```
system     = two-level atom (dim 2), Gaussian input pulse, modest optical depth
n_z=200, n_t=1000, n_v=1, fields=1, float64, n_params=50 (pulse basis)
```
Memory sweeps (n_z, n_t, dim, n_v) are run at a **memory-bound** config, not the
reference: dim 9, n_v 10. Reference cost is ~51 MB of raw rho, which won't stress anything.

| Sweep | Values | Systems for dim |
|-------|--------|-----------------|
| n_z | 50, 100, 200, 500, 1000, 2000 | |
| n_t | 200, 500, 1000, 2000, 5000 | |
| dim | 2, 4, 9, 16 | two-level, V/Lambda+1, 9-level, 16-level (define in T0.3) |
| n_v | 1, 10, 50, 100 | |

Stretch point: n_z=500, n_t=2000, dim=9, n_v=50.

## Layout

```
spike/
  DESIGN.md  PLAN.md  RESULTS.md (written at the end)
  pyproject.toml          # separate env: jax, diffrax, equinox, numpy, scipy, pytest, matplotlib
  refs/                   # ground-truth .npz from clerq (T0)
  mbjax/                  # prototype: lindblad.py, drive.py, propagate.py, strategies.py
  tests/
  bench/                  # sweep scripts, results/*.json
```

## Tasks

Legend: **Out** = artifact; **Done when** = acceptance. Phases are sequential; stop at a
phase's gate before proceeding.

### Phase 0 - Baseline and environment

- **T0.1 Environment.** `spike/pyproject.toml`, uv env, current JAX + diffrax, `jax_enable_x64`. Smoke-test `diffeqsolve` + `jax.grad` on a toy ODE. Record pinned versions. **Done when** toy grad matches analytic.
- **T0.2 Memory probe (cheap, do first).** Toy scan-of-ODE-solves; measure peak memory for plain scan vs `jax.checkpoint` body with `RecursiveCheckpointAdjoint`. **Done when** we know whether S2 ~ S3. Informs whether the strategy matrix needs to change.
- **T0.3 Define systems.** Hamiltonians, decay channels, dipole operators for dim 2/4/9/16 in a shared spec (plain dict / npz) usable by clerq and mbjax.
- **T0.4 Reference data from clerq.** Run PropagationSolver for the reference problem and a convergence study (dz, rtol halved) so "matches" is defined at converged settings. Save `refs/*.npz` with the exact z-step scheme (Euler/Adams) noted. **Done when** baseline output changes < 1e-6 under refinement.
- **T0.5 FD harness.** Generic central-FD gradient helper with fixed step, tight tolerance, step-size scan.

### Phase 1 - Inner ODE

- **T1.1** Lindblad RHS in real form for general dim from the T0.3 spec.
- **T1.2** Single-slice solve with a fixed analytic drive; compare rho(t) to clerq at z=0. **Done when** max abs diff < 1e-8.
- **T1.3** `jax.grad` w.r.t. Rabi frequency and detuning; check vs FD (T0.5).

### Phase 2 - Non-autonomous drive

- **T2.1** Drive from a sampled array: linear and cubic interpolation via `ControlTerm` (or interpolated vector field). Compare both against analytic drive.
- **T2.2** `jax.grad` w.r.t. the driving array; check vs FD and `jacfwd`. **Done when** rel. error <= 1e-5.
- **T2.3** Decide interpolation order and `SaveAt` strategy (ts vs dense).

### Phase 3 - Spatial propagation

- **T3.1** `lax.scan` over z: solve, polarisation, field update. Match z-step scheme of T0.4.
- **T3.2** Forward match vs `refs/` across the medium. **Done when** rel. error vs baseline <= 1e-6 at converged settings.
- **T3.3** `value_and_grad` of ||E_out - E_target||^2 w.r.t. a system parameter and the pulse basis; check vs FD (A5) and `jacfwd`.
- **T3.4** Velocity classes via `vmap` inside the slice body (A7). Validate vs clerq with n_v > 1.
- **Milestone:** end-to-end differentiable solver. Commit a tag.

### Phase 4 - Strategies and memory profiling (gate phase)

- **T4.1** Implement S1, S2, S3, S4, S5, FWD behind one interface (`strategies.py`).
- **T4.2** Measurement harness: subprocess per cell (clean peak RSS), forward ms, value_and_grad ms, ratio, OOM/timeout capture, JSON out.
- **T4.3** Run sweeps for each axis at the memory-bound config; run stretch point.
- **T4.4** Tuning: `jax.checkpoint` policies (`nothing_saveable`, `checkpoint_dots`, ...), inner `checkpoints=` count, chunk size for S5. Tsit5 vs Kvaerno comparison on representative systems; adaptive vs fixed step under vmap.
- **T4.5** Evaluate Gates 1-3 and the wall-clock gate; write results table.

### Phase 5 - Optimisation demo and AD-vs-FD

- **T5.1** AD vs FD vs jacfwd at n_params in {1, 10, 50, 100, 500} on reference size (Gate 4).
- **T5.2** One inverse-design problem (pulse shaping to a target output field) with L-BFGS/Adam; convergence plot; wall-clock vs FD-driven optimiser.
- **T5.3** Evaluate Gate 4.

### Phase 6 - Decision

- **T6.1** `RESULTS.md`: gate table (pass/fail with numbers), the two memory/AD-FD tables, figures, versions, caveats (CPU f64; no GPU timing).
- **T6.2** Apply the decision rule from DESIGN.md. Note which failures are fixable.
- **T6.3** If go: list the CUDA timing run and the clerq integration plan as follow-ups.

## Gates (final form)

| Gate | Criterion | No-go if |
|------|-----------|----------|
| 1 Fits | S4 or S5 completes stretch point in <= 64 GB host RAM; projected <= 32 GB device | OOM |
| 2 Scales | Peak mem(n_z=2000) <= 4x mem(n_z=200) at memory-bound config, above floor | > 10x |
| 3 Correct | AD vs fixed-step FD and vs jacfwd rel. err <= 1e-5 (f64) | > 1e-2 or unstable |
| 4 Useful | AD >= 10x faster than FD at N=100, beats FD at N=50; backward/forward <= 10x. Report vs JAX-FD **and vs clerq (PropagationSolver) FD** (reference forward: 2.3 s) | crossover above N=100 |

Decision rule unchanged: 4 pass -> build; 1 fail -> diagnose; >= 2 fail -> stop.

## Risks

- Diffrax under `vmap` + reverse mode with adaptive steps may be slow or memory-heavy: fixed-step fallback in T4.4.
- CPU timings won't transfer to GPU: reported as relative (AD/FD ratios), not absolute.
- clerq's z-stepping scheme differing from a plain scan Euler step: handled in T0.4/T3.1; may need Adams-Bashforth in mbjax.
- Stretch point may take hours on CPU: run last, in background, with a timeout.

## Status log

### T0.1 done
jax 0.11.2, diffrax 0.7.2, equinox 0.13.8, CPU, f64. `tests/test_smoke.py` passes
(grad and `jacfwd` vs analytic). `RecursiveCheckpointAdjoint` is reverse-mode only:
the FWD strategy must use `dfx.ForwardMode()` as the inner adjoint.

### T0.2 done (toy chain, Tsit5 + ConstantStepSize, n_t=500; `bench/probe_memory.py`)
Metric: XLA `memory_analysis().temp_size_in_bytes` of the compiled `value_and_grad`.
Peak RSS is not informative at these sizes (~400-600 MiB interpreter/compile floor).

| n_z | dim | direct | S2 | S3 | S5 |
|-----|-----|--------|----|----|----|
| 50  | 4 | 1.1 | 0.8 | 0.4 | 0.3 |
| 100 | 4 | 1.5 | 1.4 | 0.6 | 0.3 |
| 200 | 4 | 2.2 | 2.6 | 1.0 | 0.3 |
| 400 | 4 | 3.8 | 5.0 | 1.7 | 0.4 |
| 100 | 9 | 15.0 | 5.3 | 1.4 | 1.1 |

MiB. Findings:
- S2 and S3 **do differ**: S3 is ~2.5-4x lower, but both are linear in n_z. My expectation that
  they would be near-identical was wrong; keep both.
- S5 (chunked scan) is flat across n_z: the sqrt scaling is real and needs the nested scan.
- Everything is tiny. Extrapolating S2 linearly to the stretch point (n_t x4, n_z x5, n_v=50 if
  stored) gives roughly 5 GB: **naive S2 may already fit in 64 GB**, so Gate 1 and the
  "naive OOMs" story of Figure 1 are at risk of being trivially true/false. Extrapolation
  only; to be measured.
- Timings are single-shot and noisy; not used. Backward/forward ratio not yet measured.

Next: T0.2b, one stretch-config cell (dim 9, n_t 2000, n_z 500, vmap n_v=50) for S2/S3/S5
to settle whether the memory matrix is worth running at full scale, then T0.3.

### T0.2b done: stretch point (dim 9, n_t 2000, n_z 500, n_v 50 via vmap), toy chain
Executed cells (single run each; S5 killed after ~17 min, not measured):

| | XLA temp | peak RSS | value_and_grad | loss | grad norm |
|-|-|-|-|-|-|
| S2 | 2121 MiB | 2559 MiB | 6356 s | 416.46981527 | 39.76781311 |
| S3 | 198 MiB | 666 MiB | 7480 s | 416.46981527 | 39.76781311 |
| S5 | 192 MiB | - | - | - | - |

- Gate 1 is trivially satisfied by naive S2 on this toy (2.5 GiB at the stretch point).
- S3 uses ~10x less temp memory than S2 for ~18% more time; identical results.

### Profiling the slow stretch run (`bench/profile_slice.py`, one slice, dim 9, n_t 2000, n_v 50)

| variant | fwd (s) | grad (s) | CPU/wall (grad) |
|-|-|-|-|
| base: vmap over classes, save full state | 2.46 | 11.2 | 2.0 |
| pfn: vmap, save only polarisation | 2.45 | 11.0 | 2.0 |
| mat: state as (n, n_v) matrix, 3 GEMMs, save polarisation | 0.72 | 3.69 | 1.1 |
| mat_nosave: as mat, save final state only | 0.72 | 3.70 | 1.1 |

- 500 x base grad (11.2 s) = 5600 s, matching the observed 6356 s: the profile reproduces the run.
- Saving is NOT the cost (pfn = base, mat = mat_nosave).
- The cost is the `vmap` of the matvec: it lowers to batched small dots that burn ~2 cores for
  3x less throughput. Carrying the velocity classes as a trailing axis of the state and
  writing the RHS as plain GEMMs is ~3x faster, single-core, and gives the same forward
  value and gradient to ~1e-14 relative.
- `mat` runs at ~30 GFLOP/s, near single-core f64 peak, so it is compute-bound. Remaining
  headroom: multiple cores across velocity classes, and sparsity (see caveat).
- Caveat: the toy uses dense random 81x81 matrices. A real Lindbladian is sparse
  (commutators/anticommutators, ~O(dim^3)), so real cost per RHS should be lower.

**Amendment A7 revised:** batch Doppler classes as a trailing state axis inside the slice
body (not `jax.vmap` of the solve). Same carry (n_t), same memory behaviour. Gradient
equivalence verified on one slice (grad norm 1.47119565 both).
Expected stretch gradient with `mat`: ~3.7 s x 500 = ~31 min (S2-like), vs 1.8 h.

Next: T0.3 (define systems), T0.4 (clerq reference data).

### T0.3 done: system definitions (`mbjax/systems.py`, `systems/*.json`, `tests/test_systems.py`)
Specs are clerq atom dicts (`System(**spec)` works directly); `operators()` rebuilds H and
collapse operators in clerq's conventions without QuTiP, `liouvillian()` gives a dense
n^2 x n^2 generator for cross-checks.

| name | dim | content |
|-|-|-|
| two_level | 2 | g <-> e, one field |
| j12_12 | 4 | Jg=1/2 -> Je=1/2, x-polarised probe (q=+/-1, CG-weighted) |
| j1_2_res | 9 | Jg=1 -> Je=2 + metastable reservoir level (10% of decay) |
| j72_72 | 16 | Jg=7/2 -> Je=7/2 |

- One field per system (scalar single-field spike). Decay branching = CG^2, one collapse
  operator per (g, e); initial state is an unpolarised ground manifold.
- Correction to the plan's dim-9 idea: a pure two-manifold J system can only have even dim
  (half-integer <-> half-integer, integer <-> integer), hence the reservoir level.
- 18 tests pass: dimension, branching sums to 1 per excited level, Liouvillian is
  trace-preserving with no growing modes, and H_Omega, H_Delta, collapse operators and
  initial state all equal clerq's `System` to 1e-12.
- Note: these are the spike's *physical* systems for T0.4 reference runs. The T0.2 memory
  probes used random dense matrices; real Lindbladians are sparse, so absolute cost there
  was pessimistic.

Next: T0.4 (reference data from clerq with convergence study), T0.5 (FD harness).

### T0.4 done: clerq reference data (`refs/make_reference.py`, `refs/convergence.json`)

Reference problem (frozen): two_level, Gamma=1, Gaussian probe fwhm 1, cyclic peak Rabi
ampl=1.0 (Omega_peak = 2*pi), t in [-3, 9], square medium z in [0, 1], interaction
strength g=0.5. Transmits ~35% of input energy (peak 0.62), visibly nonlinear. Settings
n_z=200, n_t=1000, Adams-Bashforth, rtol 1e-8 -> `refs/two_level_ref.npz` (2.3 s);
converged truth (z_steps 3200, t_steps 2000, rtol 1e-10, 74 s) -> `refs/two_level_truth.npz`
(gitignored, regenerate with `uv run python refs/make_reference.py convergence`).

Convergence of the output field, relative L2 vs truth:

| z_steps | Euler | Adams-Bashforth |
|-|-|-|
| 25 | 3.5e-3 | 1.4e-4 |
| 100 | 8.8e-4 | 8.9e-6 |
| 200 | 4.4e-4 | 2.2e-6 |
| 400 | 2.2e-4 | 5.5e-7 |
| 1600 | 5.5e-5 | 2.6e-8 |

Euler is 1st order, AB 2nd order (error /4 per doubling). At z_steps=400: ODE rtol from
1e-4 to 1e-10 changes the error by <1%; t_steps 250..2000 changes it by <40%. **The error is
entirely z-discretisation.** The reference settings (n_z=200) are 2.2e-6 from the truth, so
"matches baseline" must mean *same scheme, same z_steps*, not agreement with the truth.

Conventions the JAX scan must reproduce (found by reading `propagation.py`):
1. Only Re(Omega) enters H (`Omegas.real / 2pi`); Omega itself is complex. Identical to
   the full complex coupling for symmetric detuning grids / resonant two-level, which is
   what the references use.
2. Source term: dOmega/dz = 2i * N(z_next) * (2*pi*g) * sum_k factor_k^2 * rho_{a_k b_k};
   N is evaluated at the *next* z for both AB terms.
3. AB bootstraps with a single Euler step; later steps use 1.5*f_this - 0.5*f_prev.
4. `states_zt[j+1]` is the state solved at z_j (off-by-one; states_zt[0] == states_zt[1]),
   so only `Omegas_zt` is a valid comparison target.
5. Frequencies are cyclic: collapse c = sqrt(2*pi*Gamma) sigma; H_Omega = pi*Omega_cyclic.

Decisions for the JAX implementation:
- Implement AB2 (not Euler): the carry gains one extra n_t array (previous source term),
  negligible memory. Euler would need ~4x the z_steps for the same error (2.2e-6 needs 200
  AB steps vs ~4500 Euler).
- The Phase 3 forward-match tolerance can be tight (1e-6 relative at identical z_steps),
  because the ODE tolerance and t-resolution are not the limiting error. Loose ODE rtol
  (1e-4..1e-6) suffices for forward runs; keep tight rtol only for gradient/FD checks.

Next: T0.5 (FD harness), then Phase 1.

### Note: clerq version
Developed against `origin/next` (178cb45): `MBSolve` -> `PropagationSolver`, `OBAtom` -> `System`,
`mbsolve()` -> `solve()`, `atom` -> `system`, modules `propagation.py` / `lindblad.py` /
`system.py`. The branch was first cut from a stale master (295aa89) and fast-forwarded; the
T0.4 reference output was regenerated on `next` and is bit-identical (max diff 0.0), so the
conventions above are unchanged by the rename. Old names still work with a DeprecationWarning.

### T0.5 done: FD harness and gradient correctness on the toy chain
`mbjax/fd.py` (`central_fd`, `fd_step_scan`, `check_gradient`), `tests/test_gradients.py`.
Toy chain, fixed-step Tsit5, gradient w.r.t. the input field samples through n_z=6
interpolated-drive solves (n_t=60, dim 4):

| check | result |
|-|-|
| AD vs central FD, S2/S3 (best step 1e-3..1e-4) | rel err ~1e-11 (1e-9 at h=1e-2, 1e-9 at 1e-6) |
| reverse mode (S3) vs `jacfwd` + `ForwardMode` | 2e-16 |
| direct / S2 / S3 / S5 gradients | agree to < 1e-10 |

So gradients flow correctly through the scan, the interpolated drive and the checkpointing.
Caveat: linear-in-state ODE with a quadratic loss; the real Lindblad problem (Phase 2/3)
still needs its own check.

### Adaptive-step memory (concern #2; `bench/adaptive_memory.py`, compile-time temp, MiB)

| config | strategy | const step | PID default | max_steps 4096 | 16384 | 65536 |
|-|-|-|-|-|-|-|
| ref (200/1000/dim4/nv1) | S2 | 4.2 | 4.9 | 8.0 | 14.1 | 26.3 |
| ref | S3 | 1.9 | 1.9 | 1.9 | 1.9 | 1.9 |
| stretch (500/2000/dim9/nv50) | S2 | 2120 | 2180 | 3102 | 6067 | 12028 |
| stretch | S3 | 198 | 142 | 144 | 150 | 162 |

- With adaptive stepping S2's memory is set by `max_steps` (static checkpoint buffers), so
  naive S2 only "fits in 2 GiB" if `max_steps` is small; S3 is flat at ~150 MiB regardless.
- This restores a real argument for S3: robustness to `max_steps`, not just a 10x constant.
- Caveat: static allocation, not touched pages; actual RSS may be lower (not measured).

### Open concerns after phase 0
1. Real-Lindblad gradient correctness not yet checked (Phase 2/3).
2. All timings single-run and noisy; repeat before quoting ratios.
3. Multi-level and velocity-class clerq references not yet generated (needed by T3.4).
4. Stretch-point S5 timing never measured (killed).

## Phase 1 done: inner Lindblad ODE (`mbjax/lindblad.py`, `tests/test_lindblad.py`)

Single-slice solve of the master equation in JAX/diffrax (Tsit5, PID rtol 1e-10, complex rho
carried directly), analytic Gaussian drive, detuning 0.3, t in [-2, 3].

| system | max abs diff vs clerq LindbladSolver | peak excited pop | AD vs FD (ampl, detuning) |
|-|-|-|-|
| two_level | 9e-11 | 0.25 | 1.1e-9 |
| j12_12 | 4e-11 | 0.11 | 4.8e-11 |
| j1_2_res | 1.3e-10 | 0.17 | 1.3e-10 |
| j72_72 | 3.8e-11 | 0.095 | 3.9e-10 |

- Not vacuous: a 1e-3 shift in detuning moves rho by 4e-5..6e-4, so the 1e-10 agreement
  discriminates by ~5 orders of magnitude. Trace = 1 and Hermiticity hold to 1e-8.
- Reverse mode (`RecursiveCheckpointAdjoint`) matches `jacfwd` + `ForwardMode` to < 1e-8
  and central FD to < 1e-6 (tests), actual ~1e-9..1e-11.
- Gradients are w.r.t. pulse amplitude and detuning, through a random Hermitian
  observable summed over time, so all matrix elements are exercised.
- Weak spot: the j72_72 drive is weak (rho moves only 0.035 from rho0 because the factors
  are CG/sqrt2); use a stronger drive there in Phase 3 so dim 16 is properly stressed.
- Cost: a 1000-point slice at dim 16 compiles and runs in 0.5 s.

Next: Phase 2, drive from a sampled array via interpolation, gradient w.r.t. the field.
