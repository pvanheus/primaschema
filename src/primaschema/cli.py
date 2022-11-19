import sys
import logging

import defopt

from pathlib import Path

import primaschema.lib as lib


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


def hash_bed(bed_path: Path):
    hex_digest = lib.hash_bed(bed_path)
    print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_ref(ref_path: Path):
    hex_digest = lib.hash_ref(ref_path)
    print("Reference checksum:", file=sys.stderr)
    print(hex_digest)


def validate(scheme_dir: Path):
    return lib.validate(scheme_dir)


def build(scheme_dir: Path, out_dir: Path = Path(), force: bool = False):
    lib.build(scheme_dir, out_dir=out_dir, force=force)


def main():
    defopt.run(
        {
            "hash-ref": hash_ref,
            "hash-bed": hash_bed,
            "validate": validate,
            "build": build,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
