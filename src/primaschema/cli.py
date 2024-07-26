import sys
from pathlib import Path

import defopt

from . import lib


def hash_bed(bed_path: Path):
    """
    Generate a bed file checksum

    :arg bed_path: path of bed file
    """
    hex_digest = lib.hash_bed(bed_path)
    print("BED checksum:", file=sys.stderr)
    print(hex_digest)


def hash_ref(ref_path: Path):
    """
    Generate reference sequence checksum

    :arg ref_path: path of reference sequence
    """
    hex_digest = lib.hash_ref(ref_path)
    print("Reference checksum:", file=sys.stderr)
    print(hex_digest)


def validate(scheme_dir: Path, full: bool = False):
    """
    Validate a primer scheme bundle containing info.yml, primer.bed and reference.fasta

    :arg scheme_dir: path of scheme.bed file
    :arg out_dir: path of directory in which to save primer.bed
    :arg force: overwrite existing output files
    :arg full: perform meticulous validation using full model
    """
    return lib.validate(scheme_dir, full=full)


def validate_recursive(root_dir: Path, full: bool = False, force: bool = False):
    """
    Recursively validate primer scheme bundles in the specified directory

    :arg root_dir: path in which to search for schemes
    :arg full: perform meticulous validation using full model
    :arg force: overwrite existing schemes and ignore hash check failures
    """
    lib.validate_recursive(root_dir=root_dir, full=full, force=force)


def build(
    scheme_dir: Path,
    out_dir: Path = Path("built"),
    full: bool = False,
    force: bool = False,
):
    """
    Build a primer scheme bundle containing info.yml, primer.bed and reference.fasta

    :arg scheme_dir: path of input scheme directory
    :arg out_dir: path of directory in which to save scheme
    :arg full: perform meticulous validation using full model
    :arg force: overwrite existing output files
    """
    lib.build(scheme_dir=scheme_dir, out_dir=out_dir, full=full, force=force)


def build_recursive(
    root_dir: Path, full: bool = False, force: bool = False, nested: bool = False
):
    """
    Recursively build primer scheme bundles in the specified directory

    :arg root_dir: path in which to search for schemes
    :arg full: perform meticulous validation using full model
    :arg force: overwrite existing schemes and ignore hash check failures
    :arg nested: build definitions dir structure of organism/scheme/amplicon_length/version
    """
    lib.build_recursive(root_dir=root_dir, full=full, force=force, nested=nested)


def build_manifest(root_dir: Path, out_dir: Path = Path()):
    """
    Build a complete manifest of schemes contained in the specified directory

    :arg root_dir: path in which to search for schemes
    :arg out_dir: path of directory in which to save manifest
    """
    lib.build_manifest(root_dir=root_dir, out_dir=out_dir)


def seven_to_six(bed_path: Path):
    """
    Convert a 7 column primer.bed file to a 6 column scheme.bed file by droppign a column

    :arg bed_path: path of primer.bed file
    """
    bed_str = lib.convert_primer_bed_to_scheme_bed(bed_path=bed_path)
    print(bed_str)


def six_to_seven(bed_path: Path, fasta_path: Path):
    """
    Convert a 6 column scheme.bed file to a 7 column primer.bed file using a reference sequence

    :arg bed_path: path of scheme.bed file
    :arg fasta_path: path of reference sequence
    """
    bed_str = lib.convert_scheme_bed_to_primer_bed(
        bed_path=bed_path, fasta_path=fasta_path
    )
    print(bed_str)


def diff(bed1_path: Path, bed2_path: Path, only_positions: bool = False):
    """
    Show the symmetric difference of records in two bed files

    :arg bed_path1: path of first bed file
    :arg bed_path2: path of second bed file
    :arg only_positions: Use only primer positions when computing differences
    """
    df = lib.diff(bed1_path, bed2_path, only_positions)
    if not df.empty:
        print(df.to_string(index=False))


def show_non_ref_alts(scheme_dir: Path):
    """
    Show primer records with sequences not matching the reference sequence

    :arg scheme_dir: path of input scheme directory
    """
    print(lib.show_non_ref_alts(scheme_dir=scheme_dir))


def print_intervals(bed_path: Path):
    """
    Show intervals covered by primers in a BED file

    :arg ref_path: path of bed file
    """
    all_intervals = lib.compute_intervals(bed_path)
    sorted_by_chrom = sorted(all_intervals.items())
    for chrom, intervals in sorted_by_chrom:
        sorted_interval_keys = sorted(intervals, key=lambda x: (x[0], x[1]))
        for name in sorted_interval_keys:
            interval = intervals[name]
            print(f"{chrom}\t{interval[0]}\t{interval[1]}\t{name}")


def plot(bed_path: Path, out_path: Path = Path("plot.html")):
    """
    Plot amplicon and primer coords from 7 column primer.bed
    Requires primers to be named {$scheme_id}_{$amplicon_id}_{LEFT|RIGHT}_{1|2|3…}

    :arg bed_path: path of primer.bed file
    :arg out_path: path of generated plot (with .html, .pdf, .png, or .svg extension)
    """
    lib.plot(bed_path=bed_path, out_path=out_path)


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
            "intervals": print_intervals,
            "plot": plot,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
