"""A toolkit for tiling primer scheme defintions"""

from pathlib import Path

__version__ = "0.2.0"

pkg_dir = Path(__file__).resolve().parent.parent.parent
schema_dir = pkg_dir / "schema"
primer_scheme_schema_path = schema_dir / "primer-scheme_schema.yml"
manifest_schema_path = schema_dir / "manifest_schema.json"
organisms_path = schema_dir / "organisms.yml"

# print(f"{pkg_dir=} {schema_dir=} {schema_path=}")
