import defopt

from pathlib import Path

import primaschema.lib as lib


def hash_primer_bed(bed_path: Path):
    hex_digest = lib.hash_primer_bed(bed_path)
    # print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_scheme_bed(bed_path: Path, fasta_path: Path):
    hex_digest = lib.hash_scheme_bed(bed_path, fasta_path)
    # print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_ref(ref_path: Path):
    hex_digest = lib.hash_ref(ref_path)
    # print("Reference checksum:", file=sys.stderr)
    print(hex_digest)


def build():
    pass


def validate(scheme_path: str):
    print(lib.validate_yaml(scheme_path))


def main():
    defopt.run(
        {
            "hash-ref": hash_ref,
            "hash-primer-bed": hash_primer_bed,
            "hash-scheme-bed": hash_scheme_bed,
            "build": build,
            "validate": validate,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
