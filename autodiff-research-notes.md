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
