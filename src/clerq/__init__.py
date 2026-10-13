from importlib.metadata import PackageNotFoundError, version

from clerq.field import Field
from clerq.lindblad import LindbladSolver
from clerq.propagation import PropagationSolver
from clerq.system import System

try:
    __version__ = version("clerq")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["Field", "LindbladSolver", "PropagationSolver", "System"]
