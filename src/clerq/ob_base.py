# -*- coding: utf-8 -*-

"""Deprecated module, kept so that code written for earlier versions still
imports. Use :class:`clerq.system_base.SystemBase` instead."""

from clerq import system_base as _new
from clerq._deprecation import warn_deprecated


def __getattr__(name: str):
    if name == "OBBase":
        warn_deprecated("clerq.ob_base.OBBase", "clerq.system_base.SystemBase")
        return _new.SystemBase
    try:
        return getattr(_new, name)
    except AttributeError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
