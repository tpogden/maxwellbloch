"""Spectral plot primitives: spectrum, spectrum_overlay."""

from __future__ import annotations

from typing import TYPE_CHECKING

import plotly.graph_objects as go

from maxwellbloch import spectral
from maxwellbloch.plot import theme as _theme  # noqa: F401 — registers template

if TYPE_CHECKING:
    from maxwellbloch.mb_solve import MBSolve

_TEMPLATE = "maxwellbloch"


def spectrum(
    mbs: MBSolve,
    field_idx: int = 0,
    *,
    z_idx: int = -1,
    show_dispersion: bool = False,
    width: int = 700,
    height: int = 400,
) -> go.Figure:
    """Absorption (and optionally dispersion) spectrum via Fourier transform.

    Uses the Fourier transform of the time-domain field to extract the
    frequency-domain absorption at z position *z_idx*.

    Args:
        mbs: A solved MBSolve instance. The probe field should be a weak
            pulse (so the linear-response approximation holds).
        field_idx: Index of the probe field.
        z_idx: z-step index at which to compute absorption.
        show_dispersion: If True, also plot the dispersion on a secondary axis.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    freqs = spectral.freq_list(mbs)
    absorp = spectral.absorption(mbs, field_idx, z_idx=z_idx)

    label = mbs.atom.fields[field_idx].label or f"field {field_idx}"

    fig = go.Figure(
        layout=go.Layout(
            title=f"Absorption spectrum — {label}",
            xaxis_title="frequency (γ)",
            yaxis_title="absorption (a.u.)",
            width=width,
            height=height,
            template=_TEMPLATE,
        )
    )
    fig.add_trace(go.Scatter(x=freqs, y=absorp, mode="lines", name="absorption"))

    if show_dispersion:
        disp = spectral.dispersion(mbs, field_idx, z_idx=z_idx)
        fig.add_trace(
            go.Scatter(
                x=freqs,
                y=disp,
                mode="lines",
                name="dispersion",
                yaxis="y2",
            )
        )
        fig.update_layout(
            yaxis2=dict(
                title="dispersion (a.u.)",
                overlaying="y",
                side="right",
                showgrid=False,
            )
        )

    return fig


def spectrum_overlay(
    mbs_list: list[MBSolve],
    field_idx: int = 0,
    *,
    z_idx: int = -1,
    labels: list[str] | None = None,
    width: int = 900,
    height: int = 400,
) -> go.Figure:
    """Overlay absorption spectra from multiple MBSolve results.

    Args:
        mbs_list: List of solved MBSolve instances.
        field_idx: Field index to plot in each.
        z_idx: z-step index for absorption.
        labels: Trace labels, one per element of mbs_list.
        width: Figure width in pixels.
        height: Figure height in pixels.

    Returns:
        plotly Figure.
    """
    if labels is None:
        labels = [f"run {i}" for i in range(len(mbs_list))]

    fig = go.Figure(
        layout=go.Layout(
            title="Absorption spectra",
            xaxis_title="frequency (γ)",
            yaxis_title="absorption (a.u.)",
            width=width,
            height=height,
            template=_TEMPLATE,
        )
    )
    for mbs, lbl in zip(mbs_list, labels, strict=False):
        freqs = spectral.freq_list(mbs)
        absorp = spectral.absorption(mbs, field_idx, z_idx=z_idx)
        fig.add_trace(go.Scatter(x=freqs, y=absorp, mode="lines", name=lbl))

    return fig
