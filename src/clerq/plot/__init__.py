"""clerq plotting module (Plotly-based).

Install the optional plotting dependencies with::

    pip install clerq[plot]

All primitives return a ``plotly.graph_objects.Figure``; they never call
``.show()`` internally. Use ``fig.show(renderer='notebook_connected')``
in Jupyter notebooks for interactive output, or ``fig.write_image(path)``
for static PNG export (requires kaleido).

Example::

    from clerq import PropagationSolver, plot

    mbs = PropagationSolver.from_json_str(...)
    mbs.solve()

    fig = plot.field_spacetime(mbs)
    fig.show(renderer='notebook_connected')
"""

try:
    import plotly  # noqa: F401
except ImportError as e:
    raise ImportError(
        "The clerq.plot module requires plotly. "
        "Install it with: pip install clerq[plot]"
    ) from e

from clerq.plot.fields import (
    field_envelope,
    field_profile,
    field_spacetime,
    field_z_profile_anim,
    pulse_area,
)
from clerq.plot.spectra import spectrum, spectrum_overlay
from clerq.plot.states import coherence, population, population_spacetime

__all__ = [
    "field_profile",
    "field_spacetime",
    "field_envelope",
    "field_z_profile_anim",
    "pulse_area",
    "spectrum",
    "spectrum_overlay",
    "population",
    "population_spacetime",
    "coherence",
]
