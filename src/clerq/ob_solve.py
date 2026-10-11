# -*- coding: utf-8 -*-

"""Deprecated module, kept so that code written for earlier versions still
imports. Use :class:`clerq.LindbladSolver` instead."""

from clerq import lindblad as _new
from clerq._deprecation import warn_deprecated


def __getattr__(name: str):
    if name == "OBSolve":
        warn_deprecated("clerq.ob_solve.OBSolve", "clerq.LindbladSolver")
        return _new.LindbladSolver
    try:
        return getattr(_new, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
