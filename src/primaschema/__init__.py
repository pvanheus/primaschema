from pathlib import Path

__version__ = "0.2.0"

pkg_dir = Path(__file__)
data_dir = pkg_dir.parent / "data"
schema_dir = pkg_dir.parents[2] / "schemas"
