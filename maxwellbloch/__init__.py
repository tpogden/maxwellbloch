from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("MaxwellBloch")
except PackageNotFoundError:
    __version__ = "unknown"
