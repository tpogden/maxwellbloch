"""Packaging and precision guards for the optional ``clerq.autodiff`` backend."""

import importlib
import sys

import pytest

pytestmark = pytest.mark.jax


def test_helpful_error_when_jax_is_missing(monkeypatch):
    """Importing clerq.autodiff without JAX says how to install it."""
    # None in sys.modules makes `import jax` raise ImportError, whether or not
    # JAX is installed, so this runs in every environment.
    monkeypatch.setitem(sys.modules, "jax", None)
    monkeypatch.delitem(sys.modules, "clerq.autodiff", raising=False)
    monkeypatch.delitem(sys.modules, "clerq.autodiff._config", raising=False)
    with pytest.raises(ImportError, match=r"clerq\[jax\]"):
        importlib.import_module("clerq.autodiff")


def test_import_clerq_does_not_import_autodiff():
    """The core package must not pull in the optional backend."""
    import clerq  # noqa: F401

    assert not hasattr(clerq, "autodiff") or "jax" in sys.modules


class TestFloat64Guard:
    @pytest.fixture(autouse=True)
    def _jax(self):
        self.jax = pytest.importorskip("jax")
        pytest.importorskip("diffrax")
        self.original = self.jax.config.jax_enable_x64
        yield
        self.jax.config.update("jax_enable_x64", self.original)

    def test_raises_when_x64_is_off(self):
        from clerq.autodiff import require_float64

        self.jax.config.update("jax_enable_x64", False)
        with pytest.raises(RuntimeError, match="JAX_ENABLE_X64"):
            require_float64()

    def test_does_not_change_the_global_setting(self):
        from clerq.autodiff import require_float64

        self.jax.config.update("jax_enable_x64", False)
        with pytest.raises(RuntimeError):
            require_float64()
        assert self.jax.config.jax_enable_x64 is False

    def test_passes_when_x64_is_on(self):
        from clerq.autodiff import require_float64

        self.jax.config.update("jax_enable_x64", True)
        require_float64()
