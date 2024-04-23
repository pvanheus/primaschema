import logging
import sys
from pathlib import Path

import defopt

import primaschema.lib as lib


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


def hash_bed(bed_path: Path):
    """
    Generate a bed file checksum

    :arg ref_path: Path of bed file
    """
    hex_digest = lib.hash_bed(bed_path)
    print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_ref(ref_path: Path):
    """
    Generate reference sequence checksum

    :arg ref_path: Path of reference sequence
    """
    hex_digest = lib.hash_ref(ref_path)
    print("Reference checksum:", file=sys.stderr)
    print(hex_digest)


def validate(scheme_dir: Path):
    """
    Validate a primer scheme bundle containing info.yml, primer.bed and reference.fasta

    :arg scheme_dir: Path of scheme.bed file
    :arg out_dir: Path of directory in which to save primer.bed
    :arg force: Overwrite existing output files
    """
    return lib.validate(scheme_dir)


def validate_recursive(root_dir: Path, force: bool = False):
    """
    Recursively validate primer scheme bundles in the specified directory

    :arg root_dir: Path in which to search for schemes
    :arg force: Overwrite existing schemes and ignore hash check failures
    """
    lib.validate_recursive(root_dir=root_dir, force=force)


def build(scheme_dir: Path, out_dir: Path = Path(), force: bool = False):
    """
    Build a primer scheme bundle containing info.yml, primer.bed and reference.fasta

    :arg scheme_dir: Path of input scheme directory
    :arg out_dir: Path of directory in which to save scheme
    :arg force: Overwrite existing output files
    """
    lib.build(scheme_dir=scheme_dir, out_dir=out_dir, force=force)


def build_recursive(root_dir: Path, force: bool = False, nested: bool = False):
    """
    Recursively build primer scheme bundles in the specified directory

    :arg root_dir: Path in which to search for schemes
    :arg force: Overwrite existing schemes and ignore hash check failures
    :arg nested: Build definitions inside a nested dir structure of family/version
    """
    lib.build_recursive(root_dir=root_dir, force=force, nested=nested)


def build_manifest(root_dir: Path, schema_dir: Path = Path(), out_dir: Path = Path()):
    """
    Build a complete manifest of schemes contained in the specified directory

    :arg root_dir: Path in which to search for schemes
    :arg schema_dir: Path of schema directory
    :arg out_dir: Path of directory in which to save manifest
    """
    lib.build_manifest(root_dir=root_dir, schema_dir=schema_dir, out_dir=out_dir)


def seven_to_six(bed_path: Path, out_dir: Path = Path()):
    """
    Convert a 7 column primer.bed file to a 6 column scheme.bed file by droppign a column

    :arg bed_path: Path of primer.bed file
    :arg out_dir: Path of directory in which to save primer.bed
    """
    lib.convert_primer_bed_to_scheme_bed(bed_path=bed_path, out_dir=out_dir)


def six_to_seven(bed_path: Path, fasta_path: Path, out_dir: Path = Path()):
    """
    Convert a 6 column scheme.bed file to a 7 column primer.bed file using a reference sequence

    :arg bed_path: Path of scheme.bed file
    :arg fasta_path: Path of reference sequence
    :arg out_dir: Path of directory in which to save primer.bed
    """
    lib.convert_scheme_bed_to_primer_bed(
        bed_path=bed_path, fasta_path=fasta_path, out_dir=out_dir
    )


def diff(bed1_path: Path, bed2_path: Path, only_positions: bool = False):
    """
    Show the symmetric difference of records in two bed files

    :arg bed_path1: Path of first bed file
    :arg bed_path2: Path of second bed file
    :arg only_positions: Use only primer positions when computing differences
    """
    df = lib.diff(bed1_path, bed2_path, only_positions)
    if not df.empty:
        print(df.to_string(index=False))


def show_non_ref_alts(scheme_dir: Path):
    """
    Show primer records with sequences not matching the reference sequence

    :arg scheme_dir: Path of input scheme directory
    """
    print(lib.show_non_ref_alts(scheme_dir=scheme_dir))


def main():
    defopt.run(
        {
            "hash-ref": hash_ref,
            "hash-bed": hash_bed,
            "validate": validate,
            "validate-recursive": validate_recursive,
            "build": build,
            "build-recursive": build_recursive,
            "build-manifest": build_manifest,
            "diff": diff,
            "6to7": six_to_seven,
            "7to6": seven_to_six,
            "show-non-ref-alts": show_non_ref_alts,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
