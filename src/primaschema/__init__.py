"""A toolkit for tiling primer scheme defintions"""

from pathlib import Path

__version__ = "0.2.0"

pkg_dir = Path(__file__).resolve().parent.parent.parent
schema_dir = pkg_dir / "schema"
schema_path = schema_dir / "primer_scheme.yml"
# print(f"{pkg_dir=} {schema_dir=} {schema_path=}")
