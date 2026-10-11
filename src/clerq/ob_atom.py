# -*- coding: utf-8 -*-

"""Deprecated module, kept so that code written for earlier versions still
imports. Use :class:`clerq.System` instead."""

from clerq import system as _new
from clerq._deprecation import warn_deprecated


def __getattr__(name: str):
    if name == "OBAtom":
        warn_deprecated("clerq.ob_atom.OBAtom", "clerq.System")
        return _new.System
    try:
        return getattr(_new, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
