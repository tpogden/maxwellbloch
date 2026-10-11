# -*- coding: utf-8 -*-

"""Deprecated module, kept so that code written for earlier versions still
imports. Use :class:`clerq.PropagationSolver` instead."""

from clerq import propagation as _new
from clerq._deprecation import warn_deprecated


def __getattr__(name: str):
    if name == "MBSolve":
        warn_deprecated("clerq.mb_solve.MBSolve", "clerq.PropagationSolver")
        return _new.PropagationSolver
    try:
        return getattr(_new, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
