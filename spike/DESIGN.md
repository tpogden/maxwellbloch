# Spike Design: Differentiable Maxwell-Bloch Backend

## Summary

A bounded prototype to determine whether a JAX/diffrax implementation of the coupled Maxwell-Bloch (MB) propagation equations is (a) memory-tractable on realistic problem sizes, (b) numerically correct under reverse-mode autodiff, and (c) competitively fast against finite-difference gradients at the parameter counts that matter for inverse design. The output is a go/no-go decision on building the differentiable backend.

The Clerq rename is committed independently and proceeds regardless of the spike outcome. The spike's scope is strictly the differentiable-backend question: if it passes, Clerq ships with a JAX/diffrax core and an inverse-design story; if it fails, Clerq still ships as a renamed, modernised, multi-field-capable package, and the differentiability ambition is shelved or deferred.

## Background

MaxwellBloch currently solves the coupled Maxwell-Bloch system using a NumPy/SciPy backend. The proposed Clerq backend reformulates the same system in JAX, enabling end-to-end automatic differentiation through the propagation. The motivating use case is gradient-based inverse design of atomic-vapour photonic devices: optimising input pulse shapes, medium parameters, or system Hamiltonians for target output fields.

The closest existing work is **dynamiqs** (Alice & Bob et al., JAX-based, built on diffrax), which provides differentiable Lindblad solvers for *isolated* quantum systems. dynamiqs does not address spatial propagation, which is the structural novelty here.

The headline framing — to be validated by this spike — is that autodiff yields a roughly 25× speedup at N≈100 parameters versus finite differences, despite a 3–5× per-evaluation tax, making realistic inverse design tractable for the first time.

## Problem structure

In the co-moving frame with SVEA, the MB system reduces to a sequence of ODE solves coupled through a quadrature step in z:

```
for each z-slice i:
    ρ(zᵢ, ·)  = ODESolve(Lindblad, driven by E(zᵢ, t))
    P(zᵢ, t)  = Tr[d̂ · ρ(zᵢ, t)]
    E(zᵢ₊₁,t) = E(zᵢ, t) + Δz · κ · P(zᵢ, t)
```

This is a *recurrent* computation: each slice depends on the previous one through the field. The autodiff graph is a chain of ODE solves, structurally analogous to backprop-through-time for RNNs or neural ODEs. Reverse-mode AD through the full chain yields gradients of any scalar output (e.g. ‖E_out − E_target‖²) with respect to system parameters, input pulse, or medium parameters.

There is no PDE solver. The "PDE" has been eliminated by the spatial-marching scheme; what remains is an outer loop over z and an inner ODE solve in t.

## Architecture

**Outer loop**: `jax.lax.scan` over z-slices. Required (not a Python for loop) so the entire computation traces into a single graph that reverse-mode can differentiate efficiently.

**Inner solve**: `diffrax.diffeqsolve` for the Lindblad ODE at each slice.

**Inter-slice coupling**: the ODE at each slice is non-autonomous — driven by an interpolated time-series E(zᵢ, t) that is itself an output of the previous slice. Implemented via diffrax's `ControlTerm` with a differentiable interpolation (linear or cubic; not piecewise-constant). Gradients flow through the interpolated field, so the input field at slice i−1 receives gradient contributions from the ODE solve at slice i.

**Why diffrax directly, not dynamiqs**: dynamiqs is itself built on diffrax but abstracts only the inner Lindblad solve. Wrapping it inside a custom spatial scan loop adds a layer that obscures checkpointing control and `ControlTerm` integration without saving meaningful work. Going to diffrax directly is the cleaner path for the propagation problem.

Sketch:

```python
@jax.checkpoint
def slice_step(carry, z):
    field_in, _ = carry
    rho = diffeqsolve(
        lindblad_term,
        ...,
        adjoint=RecursiveCheckpointAdjoint(),
    )
    polarisation = trace(d_op @ rho)
    field_out = field_update(field_in, polarisation, dz)
    return (field_out, ...), field_out

(_, _), field_history = jax.lax.scan(slice_step, init, z_grid)
```

## Adjoint strategy

Two coupled adjoint decisions: the inner ODE has an adjoint, and the outer scan has an adjoint. They compose.

**Inner ODE: `RecursiveCheckpointAdjoint`.** This is "discretise then optimise" — backpropagate through the solver internals using checkpointing. Diffrax's documented default and the recommended choice for nearly all problems.

**Not `BacksolveAdjoint`.** The continuous adjoint approach ("optimise then discretise") integrates a backwards-in-time ODE for the gradients. For *Hamiltonian* (Schrödinger) evolution this can work, but the Lindblad equation is dissipative — decoherence and decay drive the density matrix toward steady state, and time-reversed Lindblad dynamics are non-contractive. Numerical errors get amplified rather than damped, yielding gradients that are wrong in subtle ways. BacksolveAdjoint is a known trap for dissipative systems and stays out of the spike.

**Outer scan: `jax.checkpoint` on the scan body.** Plain `jax.lax.scan` under reverse-mode AD stores the carry state at every iteration, which is the dominant memory cost (n_z × n_t × dim²). Wrapping the scan body in `jax.checkpoint` triggers recomputation on the backward pass, reducing memory at the cost of roughly 3–4× extra forward-equivalent work.

The composition gives nested checkpointing: outer `jax.checkpoint` recomputes the whole slice (including the ODE solve) on the backward pass, and within that recomputation the inner adjoint uses its own checkpoint schedule. Expected memory scaling drops from O(n_z · n_t · dim²) to roughly O(√n_z · √n_t · dim²) under a balanced schedule.

A `ForwardMode` cross-check is worth keeping in mind for low parameter counts (≤10): forward-mode AD beats reverse-mode in that regime, and gradients computed both ways must agree.

## Implementation phases

Each phase is a self-contained validation. Don't skip ahead; the architectural risks compound.

**Phase 0 — Baseline.** Pin down a reference test problem in current MaxwellBloch (two-level atom, n_z ≈ 100, n_t ≈ 500, Gaussian pulse, modest optical depth). Save the output field as ground truth for everything that follows.

**Phase 1 — Inner ODE only.** Stand up a spike repo (separate from MaxwellBloch for now). Implement the Lindblad ODE for the two-level system in JAX/diffrax and solve a single slice driven by a fixed input field. Validate against MaxwellBloch at z=0. Take `jax.grad` w.r.t. a system parameter (Rabi frequency, detuning).

**Phase 2 — Non-autonomous driving.** Replace the fixed field with an interpolated array via `ControlTerm`. The ODE is now driven by a JAX array, so gradients can flow into it. Validate again at z=0. Take `jax.grad` w.r.t. the driving field. This is the most architecturally novel piece — most diffrax examples have analytic time-dependence.

**Phase 3 — Full spatial scan.** Wrap Phases 1–2 in `jax.lax.scan` over z. Validate the forward pass matches the Phase 0 baseline across the full medium. Take a gradient of a scalar loss w.r.t. a system parameter. End-to-end differentiable solver achieved.

**Phase 4 — Memory profiling and checkpointing.** Profile peak memory across (n_z, n_t, dim). Compare adjoint variants and `jax.checkpoint` policies on the scan body. Establish the operating regime where checkpointing buys real headroom. *This is where the gate decision is made.*

**Phase 5 — Optimisation demo.** Pick one concrete inverse-design problem (input pulse shaping for a target output, or optical depth optimisation for fidelity). Run gradient-based optimisation. Compare wall-clock against finite differences. This is the headline figure for the methods paper.

## Reference problem

One fixed configuration that all sweeps are reported relative to:

```
n_z       = 200
n_t       = 1000
dim       = 4    (two-level + a small set of Doppler classes)
n_v       = 1    (single velocity class for the reference; sweep separately)
n_params  = 50
fields    = 1
precision = float64 / complex128
```

Small enough to run dozens of configurations in an afternoon; large enough that memory effects are visible; physically meaningful (close to a Phase A magneto-optical rotation or simple EIT problem). One stretch point — approximating B4 simultons (n_z=500, n_t=2000, dim=9, n_v=50) — to confirm the projections from 1D sweeps.

## Experiments

**Sweep axes** (one at a time, others held at reference):

- n_z ∈ {50, 100, 200, 500, 1000, 2000} — dominant memory axis
- n_t ∈ {200, 500, 1000, 2000, 5000} — inner ODE trajectory length
- dim ∈ {2, 4, 9, 16} — Hilbert-space scaling, capturing two-level → multi-level Λ → Rydberg-scale
- n_v ∈ {1, 10, 50, 100} — Doppler velocity classes (relevant for thermal vapour problems)

**Strategies** (in order of increasing engineering complexity):

| Strategy | Description |
|----------|-------------|
| S0 | NumPy/SciPy baseline (current MaxwellBloch, forward only) |
| S1 | JAX forward only, no AD (pure JIT cost vs S0) |
| S2 | Naive AD: `RecursiveCheckpointAdjoint` inner + plain scan |
| S3 | Scan-checkpointed: S2 + `jax.checkpoint` on scan body |
| S4 | Doubly-checkpointed: S3 + tuned inner adjoint checkpoints |

**Metrics per cell**:

- Peak memory (MiB)
- Forward wall-clock (ms)
- Backward (or `value_and_grad`) wall-clock (ms)
- Backward/forward ratio
- Gradient relative error vs central finite differences (on a problem small enough to make FD feasible)

**Two output tables**:

1. *Memory matrix* (one per sweep axis): rows = axis values, columns = strategies, cells = peak memory. OOM cells flag where each strategy falls over.
2. *AD-vs-FD scaling* at fixed reference size: rows = n_params ∈ {1, 10, 50, 100, 500}, columns = FD time / AD time / speedup. The crossover point validates or invalidates the headline framing.

## Gate criteria

Committed in advance, before any data is collected.

**Gate 1 — Fits.** Doubly-checkpointed strategy (S4) fits the stretch point in ≤32 GB of GPU memory. *No-go if S4 OOMs.*

**Gate 2 — Scales.** Memory under S3 grows as O(n_zᵅ) with α ≤ 0.7 across the n_z sweep. *No-go if α ≥ 1.0 — checkpointing isn't actually helping.*

**Gate 3 — Correct.** Autodiff gradients match central finite differences with relative error ≤1e-5 on the reference problem (f64 throughout). *No-go if gradients are unstable or wrong by more than ~1%.*

**Gate 4 — Useful.** AD beats FD at n_params=50 with ≥10× speedup at n_params=100. *No-go if the crossover is above n_params=100 — inverse design has no efficiency story.*

**Decision rule**:

- All four pass → build the differentiable backend in Clerq; draft the methods paper around the memory and speedup figures.
- One fails → diagnose. Some failures are fixable (Gate 1 may be solved by reparameterisation), some are fundamental (Gate 4 failure kills the headline claim).
- Two or more fail → strong stop signal on the differentiable backend. Clerq still ships without it; the methods paper is shelved or deferred. Reconsider scope.

## Outputs

The same spike data, recombined, becomes the empirical core of the methods paper:

- **Figure 1**: Memory scaling with n_z, four lines (S1/S2/S3/S4), log-log axes. Story: naive AD scales linearly and OOMs by n_z=500; doubly-checkpointed scales sub-linearly and remains feasible at n_z=2000.
- **Figure 2**: AD vs FD wall-clock at fixed problem size, sweeping n_params. Story: AD pays a constant overhead but wins decisively above N≈10; at N=100 it is ≈25× faster, making realistic inverse design tractable.

The gate-decision artifact *is* the paper's empirical core. No additional experimental work is required to write the methods paper if the gates pass.

## Out of scope

Explicitly not in this spike:

- Multi-field (>2) propagation. Scalar single-field validation only; multi-field structural extension comes after the gate decision.
- Counter-propagation (Phase A0 of the MaxwellBloch expansion plan). Co-propagating only.
- Phase C extensions (density-dependent mean-field, z-dependent detuning, quantum noise, DC B-field/Zeeman, Plotly plotting).
- Migration of any existing MaxwellBloch examples. Spike runs in a separate repo.
- API design or naming. The spike answers "is this feasible," not "what does the public interface look like."
- GPU optimisation beyond the default JAX behaviour.

## Open questions

To be resolved during the spike:

- Optimal `jax.checkpoint` policy: `checkpoint_dots`, `nothing_saveable`, `dots_with_no_batch_dims_saveable`, or custom. Treated as a tuning step inside S4 rather than a separate strategy.
- Whether Doppler-class batching should be inside the scan carry or as an outer `vmap`. Affects memory layout under reverse-mode.
- Whether `Tsit5` (recommended general-purpose) or `Kvaerno{3,4,5}` (stiff-friendly) is the right default solver for representative Lindblad systems. Likely Tsit5 at f64 with adaptive PIDController, but worth a small comparison.
- Whether `ForwardMode` cross-checks should be standardised at low n_params as a permanent regression test.
- Behaviour under JAX 0.5+ with the most recent diffrax; the API has been stable but worth a smoke test before committing to versions.

## Pre-spike reading

A few hours each, all worth doing before writing code:

- dynamiqs source: how they handle Lindblad + diffrax + adjoint composition.
- diffrax adjoint docs: `RecursiveCheckpointAdjoint`, `BacksolveAdjoint`, `ImplicitAdjoint`, `ForwardMode` trade-offs.
- Back-of-envelope memory estimate for the reference problem under naive S2 versus checkpointed S3.

## Status

Design draft. Not yet started.
