# Autodiff research notes (working draft, uncommitted)

Date: 2026-10-19. Branch: claude/autodiff-research-plan-449fd8. Companion to the plan in
`~/.claude/plans/users-tpo-downloads-diffable-maxwellblo-merry-charm.md`.

## 1. Literature check (4 web searches, not a systematic review)

No differentiable Maxwell-Bloch / atomic-vapour / EIT-memory tool found. Related work:
- Memory optimal control without AD: Gorshkov et al. time-reversal iteration; Phillips, Gorshkov,
  Novikova (arXiv:0805.3348); Novikova et al. (arXiv:0805.1927). Beil et al., PRA 2010 used an
  evolutionary search for EIT storage pulses.
- Differentiable quantum dynamics (single point, no propagation): dynamiqs, qutip-jax, qiskit-dynamics;
  memory-efficient GRAPE in JAX (arXiv:2203.12717, 2210.08378).
- Photonics inverse design with AD/adjoint FDTD (e.g. time-reversal differentiation of FDTD,
  Graz; Hybrid AD + adjoint, Mach. Learn.: Sci. Technol. 2024): the closest methodological
  relatives, worth citing for the "adjoint through a time-marching solver" framing.
Absence of evidence only. Before claiming "first", do a proper Scholar/arXiv pass and ask people.

## 2. Memory prototype (scratchpad `memproto.py`, not in repo)

Setup: `mbs-lambda-eit-slow-light` store configuration, interaction strength 5e2, gaussian density,
z_steps 70, z_steps_inner 20, plus excited-state decay (rate 1, both channels) and ground dephasing
(`[2,2]`, rate 0.01). Loss = probe energy out / in. Control = 141 free samples of the coupling field
via `omega0`, passed through 5*sigmoid so the peak cannot exceed the baseline (5). Adam, lr 0.1, 150 steps.

- Baseline (ramp off/on) efficiency 0.338. Optimised 0.379 (1.12x). From the clipped sigmoid init (0.301) it is 1.26x.
- That is under the 1.3x Gate B bar, with a fixed storage delay and a cap equal to the baseline peak. Not a headline yet.
- Cost: about 9 s per gradient step, 1360 s total on CPU.

Caveats found:
- The z-march is unstable at coarse settings (z_steps 40, inner 4 gave energy ratios ~1e9). Inherited from
  clerq's scheme; examples need a stability guard or a documented minimum `z_steps_inner`.
- The finite-difference spot check I ran did NOT agree (AD -2.25e-3 vs FD -1.3e-4 at one index), but the
  FD used the adaptive solver, which is noise-limited. The repo's own tests use fixed `n_steps`. So this
  check was inconclusive. UPDATE: redone with fixed-step FD (n_steps=2003, N=50, 100-sample coupling control): AD
  matches FD to ~1e-8 relative at the four largest components, so the gradients are right.
  The sweep (N=50/150/500 x delay 1/2/4 x free probe) is running overnight; first free-probe run gave efficiency
  3.89 (>1), treated as an artifact until the waveform is inspected.
- Loss was only 0.04 above its start by step 50 onwards; unconverged vs. converged landscape unknown.
- Not yet explored: longer/free delay, input-pulse shaping, higher optical depth, power-budget Pareto front.

## 3. Next
1. Redo the gradient check with fixed `n_steps` FD (and compare forward vs reverse).
2. Sweep delay and optical depth; add input-pulse parameters; compare with Gorshkov-optimal baseline.
3. If gain stays near 1.1-1.3x, lead with the spectrum-fit and sensitivity examples and keep memory as a demonstration.

## 4. Candidate datasets for the spectrum-fit example (Gate C input)

Found by search; contents not yet inspected, nothing downloaded.
- **Higgins & Hughes, Durham** — "EIT in a V-system with 87Rb vapor in the hyperfine Paschen-Back regime"
  [dataset], doi:10.15128/r2t148fh16h, CC BY 4.0, a ~220 kB zip (contents not listed on the record page). The
  record compares against a time-dependent Lindblad model averaged over interaction times, for a range of pump
  Rabi frequencies, which is close to clerq's model. Needs: clerq to handle 0.6 T Zeeman/PB level structure
  (hyperfine module + autodiff unverified), and the companion paper to be identified.
- **NIST PDR** Rydberg EIT/Autler-Townes sets: mds2-3102 (Cs, Zeeman-resolved AT splitting; "data behind every
  figure"), mds2-4115 (fundamental linewidth limit of thermal Rydberg ladder EIT), mds2-3005 (dual-tone RF
  dressed EIT, CSV). Licences not checked. Natural for the ladder/Rydberg notebook.
- ElecSus test data (github.com/jameskeaveney/ElecSusTestData): weak-probe D-line spectra; linear-response only,
  so it would test the fit machinery but not the propagation physics.
Next: ask before downloading the Durham zip (and check the other records' licences).

## 5. Reading list

Memory optimal control: Gorshkov et al. PRA 76 033804 (2007) and companion papers; Phillips, Gorshkov, Novikova
(arXiv:0805.3348); Novikova et al. (arXiv:0805.1927); Beil et al. PRA 2010 (learning-loop pulse search).
Differentiable quantum dynamics: dynamiqs, qutip-jax, qiskit-dynamics; GRAPE memory, Lu, Joshi, Dinh, Koch
(arXiv:2304.06200; J. Phys. Commun. 2024): argues AD-GRAPE is memory-hungry, so our checkpointing results
should be positioned against it; Argonne memory-efficient GRAPE (arXiv:2203.12717, 2210.08378).
Photonics adjoint/AD: time-reversal FDTD differentiation (Graz); hybrid AD+adjoint (Mach. Learn.: Sci. Technol.
2024); Colburn et al. Commun. Phys. 2021.
Thermal-vapour models: Higgins/Hughes EIT modelling; arXiv:2111.13408 (Rydberg EIT spectrum in a thermal vapour),
arXiv:2307.08545 (hyperfine Paschen-Back EIT).
Not found by search (still to locate by author name): Hush/Hedges spin-wave memory control papers.

## 6. Design: per-channel decay rates as differentiable parameters

State today: only the global `decay_scale` multiplies every rate. `build_operators` bakes the collapse operators into
`SparseDissipator(a, b, w, r)`, with `w_k = |f_k|^2` (rate and factor included), `r = bincount(b, w)`
(`_operators.py`), and zero-rate operators are filtered out.
Key observation: the dissipator is linear in `w` (`D(rho) = sum_k w_k D_k(rho)`), so making each rate a parameter is
cheap and the gradient is well defined at any rate, including zero (no sqrt singularity).
Proposal:
1. New parameter `decay_rates`, shape `(n_decays,)`, cyclic units, default from `system.decays[i]["rate"]`, one entry per
   decay dict (all its channels share the rate; `factors` stay fixed). This matches how users write JSON.
2. Keep an index `decay_of_op[k]` mapping each collapse operator to its decay dict; build `w = 2 pi rate[decay_of_op] * f_k^2`
   in JAX and `r` via `segment_sum`, replacing the baked numpy `w`, `r`.
3. Do not drop zero-rate channels at extraction (so a rate can be optimised from 0); drop only structurally empty ops.
4. Dense fallback: scale each unit-normalised operator by `sqrt(w_k)` inside the dissipator, or use the linear form directly.
5. `decay_scale` stays as a global multiplier; rates enter the cache key/shape checks like other parameters.
6. Tests: match `solve()` after changing a rate (including a `[a, a]` dephasing channel and ground-state decay);
   reverse = forward = fixed-step FD for the rate parameters; a zero-rate gradient test; bench the cost (should be ~0).
Open questions: should per-channel rather than per-dict rates be exposed; units/validation for negative rates (clip or
raise); how this interacts with `decay_scale` in `default_params`.
Effort estimate: small (days), the most valuable feature gap for the fit and memory examples.

## 7. Memory prototype: earlier results INVALIDATED (2026-10-20)

The table previously here (0.338 -> 0.364 / 0.580 at N=500; 0.493 / 0.762 at N=150) is not a memory result.
The metric counted any probe energy after the coupling returns (t >= t_on), and the optimiser exploited that:
- N=500: the optimised coupling stayed ~5 throughout, including the hold window, so nothing was stored; it was slow-light transmission.
- N=150: 87% of the optimised input probe energy arrived after t=4 and 67% after t=6, i.e. the input was simply transmitted late.
Likewise the earlier sweep's 107x "efficiencies" were near-Nyquist artifacts (section 2). What stands: the gradient check against
fixed-step FD (~1e-8), and the finding that unconstrained free-sample controls exploit the time grid.
Lesson for the paper/examples: a memory objective must fix the protocol (coupling forced to zero during the hold, input probe
supported before the hold) and count only retrieval after the coupling returns; report leakage separately.

Reference bound: Gorshkov et al. kernel K(z,z') = (d/2) exp(-d(z+z')/2) I0(d sqrt(zz')) on [0,1]; eta_s(d) = largest eigenvalue
(checked: -> 1 as d grows, (1-eta_s)*d -> ~2.9). eta_s = 0.81 (d=10), 0.92 (d=30), 0.97 (d=100). Total storage+retrieval
efficiency <= eta_s(d). clerq's optical depth is d = 2 g (g = interaction_strengths; measured 1.958 g with a narrowband probe).
Corrected protocol script: scratchpad mem3.py (uniform medium, forced hold, Doppler optional).

## 8. Corrected memory results (forced hold; 2026-10-21)

Protocol: uniform medium on [0,1], coupling forced to ~0 during the hold (t=3..5), input probe supported before the hold,
only probe energy at t >= 5 counts as retrieved, leakage reported separately. Coupling capped at 5, ground dephasing 0.01,
12 coupling + 10 probe Gaussian basis functions, 80 Adam steps (lr 0.1). Robust to finer z (and, without Doppler, finer t) within 0.3%.

| case | d | Gorshkov eta_s(d) | retrieved | leak | notes |
|-|-|-|-|-|-|
| no Doppler | 40 | ~0.94 | 0.525 | 0.071 | |
| no Doppler | 100 | 0.975 | 0.619 | 0.034 | |
| Doppler, 9 classes over +-4 widths (width 2) | 40 (nominal g) | n/a | 0.100 | 0.103 | 35 min for 80 steps (~26 s/step) |

- The Ramp/Gaussian baseline stores ~nothing in this setup (retrieved ~0, 50-90% of the pulse leaks before the hold), so no
  gain-over-baseline ratio is meaningful; the bound is the right reference. A tuned hand protocol would be a better baseline.
- Gap to bound is unexplained: candidates are unconverged optimisation, coupling cap 5, finite bandwidth (bound is adiabatic),
  dephasing 0.01, basis size. Total storage+retrieval is <= eta_s; the forward-retrieval eta_s^2 relation was not verified.
- Doppler result is NOT interpretable yet: (a) I have not verified that clerq's sign conventions preserve two-photon
  resonance for co-propagating Lambda fields; (b) the nominal d=2g no longer equals the resonant OD once classes are weighted;
  (c) the loss is not converged (last training value 0.103 vs 0.100 after the final update; started ~1e-6).
Next checks: Doppler-broadened weak-probe EIT transmission vs no-Doppler (two-photon cancellation); measure effective OD with
the classes; longer/scheduled optimisation; Gate B decision needs a tuned baseline and a converged Doppler number.

## 9. Doppler checks and longer runs (2026-10-21)

Checks (weak narrowband probe, fwhm 12, nominal d=40):
- Two-photon Doppler cancellation holds with clerq's default field settings: intensity transmission with coupling 5 is 0.968 (no Doppler) and
  0.549 (Doppler, 11 classes); flipping the coupling's `factor_doppler_shift` to -1 or its `detuning_positive` drops it to 0.007.
- Effective optical depth with the classes: d_eff = 1.246 g (vs 1.958 g without Doppler). `thermal_delta_steps=9` gives 11 classes
  (10 on linspace plus a zero class), spacing 0.89 (cyclic), finer than the linewidth. So the g=20 Doppler case has d_eff ~ 25.

Longer runs (200 Adam steps, cosine lr 0.1 -> 0.005, forced hold t=3..5, same basis, one initialisation each):
| case | d | retrieved | leak | notes |
|-|-|-|-|-|
| no Doppler | 25 | 0.460 (plateau 0.4595) | 0.102 | 23 min; z/t refinement within 0.4% |
| Doppler, 11 classes (width 2) | d_eff ~ 25 | 0.118 (still creeping up, 0.1166 -> 0.1175 over the last 20 steps) | 0.073 | 86 min (~26 s/step); z-refined 0.118 |

- Gorshkov eta_s(25) is roughly 0.90, so no-Doppler reaches ~51% of the bound; the Doppler run reaches ~26% of the no-Doppler value.
  The bound is for homogeneous broadening, so it is only a loose upper bound for the Doppler case.
- Single initialisation per case; local optima not ruled out. No tuned hand baseline exists yet, so Gate B (>=1.3x over standard protocols
  under realistic dephasing) can't be judged. The default Gaussian + ramp protocol stores nothing here.
- What is shown: the AD pipeline optimises a 22-parameter memory protocol through an 11-class Doppler medium; the physics result so far is
  only that Doppler broadening costs a large factor in this parametrisation.
Next: tuned baseline by brute-force scan over a few parameters (pulse width, timing, coupling amplitude), multiple restarts, a larger
coupling cap / basis, and a timing comparison against finite differences for the same problem.

## 10. Tuned baseline, restarts, larger cap, FD timing (2026-10-22)

Setup as section 9 (d=25 no Doppler; Doppler 11 classes with d_eff~25; 22 parameters = 12 coupling + 10 probe basis coefficients; 120 Adam steps
cosine lr). Baseline = best of 100 Gaussian probes (fwhm in 0.8..2.5, centre -1..1, constant coupling amplitude 2..10 outside the hold).

FD timing for the 22-parameter gradient (single machine, CPU, best-effort single measurements):
| problem | forward | reverse grad | grad/forward | forward-difference grad (23 fwd) | FD / reverse |
|-|-|-|-|-|-|
| no Doppler | 0.82 s | 4.7 s | 5.7x | 18.8 s | 4.0x |
| Doppler (11 classes) | 3.2 s | 20.5 s | 6.4x | 75.1 s | 3.7x |
AD and FD gradients agree to 2e-4 and 3e-5 relative (adaptive solver, so noise-limited; not a validation). Consistent with the spike: cost ratio ~ P/6.
Implication for runs: 200 Doppler steps took 86 min with AD; the same with FD would be ~5 h (estimated from the measured ratio, not run).

No Doppler, d=25 (retrieved energy fraction):
| start | cap 5 | cap 10 |
|-|-|-|
| best Gaussian scan (w=2.5, c=1, A=2: 0.268) then optimised | 0.272 | 0.274 |
| random init, seed 1 | 0.416 | 0.425 |
| random init, seed 2 | 0.377 | - |
| earlier default init, 200 steps (section 9) | 0.460 | - |
- The scan optimum sits at three grid edges (w, c, and the lowest amplitude), so the "tuned baseline" is not converged; extend the grid before quoting a ratio.
- Large spread across starts (0.27 to 0.46): local optima / slow convergence are real. Raising the cap 5 -> 10 changed little (+0.01).
- Indicative: free-form optimised 0.38-0.46 vs best-scanned Gaussian 0.27 (about 1.4-1.7x), subject to the edge caveat.
Doppler (d_eff~25), started from the scan optimum (0.0015): 0.016 (cap 5) and 0.023 (cap 10) after 120 steps, still climbing, with large leakage (0.38-0.57).
The earlier default-init 200-step run reached 0.118 (section 9), so these runs are unconverged and the start point matters a lot; no Doppler conclusion yet.
Not yet done: wider scan grid, longer/multi-start Doppler optimisation, any comparison with the analytic optimal-control protocol shape.

## 11. Tuned baseline settled (2026-10-22), Gate B assessment

Tuned baseline = Gaussian probe (fwhm w, centre c, then truncated by the input-window taper at t=2.4 and renormalised) + constant coupling amplitude A
outside the hold; 3 parameters; grid then Nelder-Mead. Note: the best centre lies beyond the window (c~3), i.e. the input is a rising exponential-like
leading edge, the known good storage shape. My first grid stopped at c=1.5 and missed it (so the earlier "1.4-1.7x" over the scan was an artifact of
a bad baseline).

| case (d~25) | tuned 3-param baseline | free-form optimiser (22 params) | ratio |
|-|-|-|-|
| no Doppler | 0.469 (grid 504 evals + NM 150; w=0.89, c=3.09, A=3.82) | 0.460 best (200 steps); 0.38-0.43 other starts | ~1.0 |
| Doppler, 11 classes | 0.067 (NM, 60 evals/start, 3 starts: 0.067, 0.064, 0.049) | 0.118 (200 steps, one start) | ~1.8 |

Gate B (>=1.3x over a standard protocol under realistic dephasing): NOT met without Doppler; tentatively met with Doppler, but the Doppler baseline
is unconverged (60 evals), the optimiser has one run, and the free-form basis (Gaussians of width ~0.67) cannot represent the sharp input cut-off, so
the no-Doppler tie may reflect the basis rather than the method. Interpretation: gradient-based free-form design earns its cost where simple
families fail (inhomogeneous broadening), not where a 3-parameter scan suffices; at 3 parameters Nelder-Mead (3 min) beat 22-parameter Adam (20 min).
To confirm before any claim: longer NM and a richer baseline family on Doppler; several optimiser starts and a finer basis; then compare on
a second optical depth.

## 12. Doppler advantage: confirmation runs (2026-10-23)

Doppler, 11 classes, d_eff~25, coupling cap 5 for everything, retrieved energy fraction (forced hold, same metric as above).
- Richer baseline: 5-parameter family (Gaussian-like probe with width, centre, shape exponent; separate write/read coupling amplitudes), Nelder-Mead 300 evals
  from two starts: 0.0681 and 0.0681 (converged to the same optimum; w=0.29, c=2.90, A_write=4.4, A_read=5 (clipped), exponent 1.52). ~17 min each.
- Free-form optimiser, finer basis (20 coupling + 24 probe coefficients = 44 parameters), 100 Adam steps, ~39 min each:
  random start 1: 0.1318; random start 2: 0.1320 (agree to 0.2%; still creeping up: 0.1301 -> 0.1318 over the last 20 steps); leak ~0.044.
  Earlier coarse-basis (22 parameters, 200 steps): 0.118.
- Ratio: ~1.9x over the best 5-parameter baseline. Gate B (>=1.3x) is met for Doppler under these conditions; without Doppler the ratio is ~1.0.
- Caveat/anomaly: a warm start from the baseline optimum ended at 0.027 (below the baseline). The basis could not reproduce the sharp input cut-off
  (initial fit 0.049 vs 0.068), and Adam at lr 0.05 first fell to 0.009 and recovered slowly. So the landscape has separate basins and the optimiser result
  depends on its start; the two random starts agreeing is the better evidence.
- Not covered: a single optical depth and Doppler width; baseline families beyond the 5-parameter one (e.g. shaped coupling edges) could narrow the gap;
  no theoretical bound for inhomogeneous broadening; waveforms not yet inspected for physical interpretation; FD optimisation not run
  (estimate at 44 parameters: forward-difference gradient ~45 forwards vs ~6 for reverse, i.e. ~7x, extrapolated from the 22-parameter timing).
Current story for a methods/application paper: gradient-based design earns its cost in the Doppler-broadened regime where low-parameter protocols are
~2x worse, and not where a 3-parameter scan suffices.

## 13. Spectrum-fit scoping: detuned-probe scans are not valid in clerq (2026-10-23)

Prototype probe (scratchpad fit0.py/fit1.py): `jax.vmap` over the probe detuning and `jax.jacrev` of per-detuning transmission both work through
`AutodiffPropagation.propagate` (25-point spectra and their Jacobians are feasible: ~12 s for a 3-point 3-parameter Jacobian on a small problem).
But a two-level check (coupling off, g=10, weak narrowband probe) gives, versus the expected Lorentzian transmission exp(-OD/(1+4 d^2)):
d=0.25: 0.062 (expected 1e-7); d=1: 0.98 (0.018); d=3: 3.7 (0.58). Transmission above 1 is unphysical.
Cause (consistent with the code): `propagation.py:670` drives the Hamiltonian with `Omegas_z[f_i].real`, i.e. only Re(Omega), which is exact on resonance
but mixes +d and -d for a detuned field. clerq's spectral-analysis notebook avoids this by Fourier-analysing a broadband resonant pulse.
Consequences:
- A spectrum-fit example must either use that FFT route (valid only for spectra symmetric about the probe centre, so no two-photon offsets, dispersion
  asymmetry or Paschen-Back multiplets such as the Durham dataset) or clerq needs a complex-envelope drive (an RWA-correct H with Omega and Omega*),
  in both the QuTiP path and `autodiff/_operators.py`/`_lindblad.py` (which integrate a real (A, B) pair).
- OPEN RISK for the Doppler memory results: velocity classes are detuned fields too. Symmetric class grids may make the Re(Omega) approximation exact on
  average (per the earlier code survey), but this has not been verified. Check against a complex-drive reference before quoting the Doppler numbers.
