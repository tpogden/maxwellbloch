"""MaxwellBloch plotting module (Plotly-based).

Install the optional plotting dependencies with::

    pip install maxwellbloch[plot]

All primitives return a ``plotly.graph_objects.Figure``; they never call
``.show()`` internally. Use ``fig.show(renderer='notebook_connected')``
in Jupyter notebooks for interactive output, or ``fig.write_image(path)``
for static PNG export (requires kaleido).

Example::

    from maxwellbloch import mb_solve, plot

    mbs = mb_solve.MBSolve.from_json_str(...)
    mbs.mbsolve()

    fig = plot.field_spacetime(mbs)
    fig.show(renderer='notebook_connected')
"""

try:
    import plotly  # noqa: F401
except ImportError as e:
    raise ImportError(
        "The maxwellbloch.plot module requires plotly. "
        "Install it with: pip install maxwellbloch[plot]"
    ) from e

from maxwellbloch.plot.fields import field_envelope, field_spacetime, pulse_area
from maxwellbloch.plot.spectra import spectrum, spectrum_overlay
from maxwellbloch.plot.states import coherence, population, population_spacetime

__all__ = [
    "field_spacetime",
    "field_envelope",
    "pulse_area",
    "spectrum",
    "spectrum_overlay",
    "population",
    "population_spacetime",
    "coherence",
]
