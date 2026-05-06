####################
  Two-level
####################

*************************************
Linear Absorption & Beer-Lambert Law
*************************************

Weak probe pulses and continuous-wave fields in a two-level absorber.
Demonstrates Beer-Lambert attenuation, the role of spontaneous decay, and
how increasing optical depth reshapes the transmitted pulse.

.. toctree::
   :maxdepth: 1

   mbs-linear-absorption

***************************************************
Self-Induced Transparency & Optical Solitons
***************************************************

Strong resonant pulses whose area exceeds π can form area-conserving solitons
via self-induced transparency (McCall–Hahn). 2π solitons propagate without
loss; larger areas break up into multiple solitons. Odd-π inputs are unstable
fixed points of the area theorem.

.. toctree::
   :maxdepth: 1

   mbs-sit-area-theorem

***************************************************
Photon Echo & Coherence Time ($T_2$) Measurement
***************************************************

A π/2 pulse followed by a π pulse creates a photon echo at $t = +\tau$ in an
inhomogeneously broadened (Doppler-broadened) medium. The echo amplitude decays
as $e^{-2\tau/T_2}$, allowing the **homogeneous** coherence time $T_2$ to be
measured even when the lineshape is dominated by inhomogeneous broadening.

.. toctree::
   :maxdepth: 1

   mbs-two-photon-echo

**************************************************
Stimulated Echo & Population Lifetime ($T_1$)
**************************************************

Three π/2 pulses with delays τ (between pulses 1–2) and $T_w$ (between
pulses 2–3) generate a stimulated photon echo at $t = +\tau$ after the
third pulse. The echo amplitude decays as
$e^{-2\tau/T_2}\,e^{-T_w/T_1}$: a $T_w$ sweep isolates $T_1$,
complementing the $T_2$ measurement from the two-pulse echo. The defining
feature is that coherence information is stored as a **population grating**
during $T_w$, so no macroscopic polarisation exists in that interval.

.. toctree::
   :maxdepth: 1

   mbs-two-stimulated-echo

