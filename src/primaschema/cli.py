import sys
import defopt

from pathlib import Path

import primaschema.lib as lib


def hash_bed(bed_path: Path):
    hex_digest = lib.hash_bed(bed_path)
    print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_primer_bed(bed_path: Path):
    hex_digest = lib.hash_primer_bed(bed_path)
    print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_ref(ref_path: Path):
    hex_digest = lib.hash_ref(ref_path)
    print("Reference checksum:", file=sys.stderr)
    print(hex_digest)


def build():
    pass


def main():
    defopt.run(
        {
            "hash-bed": hash_bed,
            "hash-ref": hash_ref,
            "build": build,
            "hash-primer-bed": hash_primer_bed,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
