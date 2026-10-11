# -*- coding: utf-8 -*-

"""Helpers for names kept only for backwards compatibility."""

import functools
import warnings
from collections.abc import Callable


def warn_deprecated(old: str, new: str, stacklevel: int = 3) -> None:
    """Warn that ``old`` is deprecated in favour of ``new``.

    Args:
        old: The deprecated name.
        new: The name to use instead.
        stacklevel: Frames up to the user code the warning should point at,
            counted from the caller of this function (the default suits a
            direct call from the deprecated function or property).
    """
    warnings.warn(
        f"{old} is deprecated and will be removed in a future release; "
        f"use {new} instead.",
        DeprecationWarning,
        stacklevel=stacklevel,
    )


def resolve_system_arg(
    system: dict | None, atom: dict | None, stacklevel: int = 4
) -> dict | None:
    """Return the ``system`` argument, accepting the deprecated ``atom``.

    Raises:
        TypeError: If both ``system`` and ``atom`` are given.
    """
    if atom is None:
        return system
    if system is not None:
        raise TypeError("Pass 'system' only; 'atom' is its deprecated name.")
    warn_deprecated("The 'atom' argument", "'system'", stacklevel=stacklevel)
    return atom


def rename_atom_key(json_dict: dict, stacklevel: int = 4) -> dict:
    """Rename the deprecated top-level ``atom`` key of a problem definition."""
    if "atom" in json_dict:
        if "system" in json_dict:
            raise ValueError(
                "Problem definition has both 'system' and 'atom' keys; "
                "'atom' is the deprecated name for 'system'."
            )
        warn_deprecated("The 'atom' key", "'system'", stacklevel=stacklevel)
        json_dict = dict(json_dict)
        json_dict["system"] = json_dict.pop("atom")
    return json_dict


def renamed_arg(old: str, new: str) -> Callable:
    """Decorator accepting the deprecated keyword ``old`` for argument ``new``."""

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if old in kwargs:
                if new in kwargs:
                    raise TypeError(
                        f"Pass '{new}' only; '{old}' is its deprecated name."
                    )
                warn_deprecated(f"The '{old}' argument", f"'{new}'")
                kwargs[new] = kwargs.pop(old)
            return func(*args, **kwargs)

        return wrapper

    return decorator
