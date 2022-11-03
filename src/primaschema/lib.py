import hashlib
from pathlib import Path

import defopt
import pandas as pd
from Bio import SeqIO


def hash_bed(bed_path: Path) -> str:
    """Legacy coordinate hashing function"""
    scheme_bed_fields = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "poolName",
        "strand",
    ]
    fields_to_hash = ["chrom", "chromStart", "chromEnd", "strand"]

    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=scheme_bed_fields,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            strand=str,
        ),
    )

    # Normalise and hash bed records
    df_norm = df.applymap(
        lambda x: x.strip() if isinstance(x, str) else x
    )  # Strip trailing and leading whitespace
    df_norm = df_norm[[*fields_to_hash]].sort_values("chromStart")  # Sort by start pos
    df_norm_text = df_norm.to_csv(sep="\t", header=False, index=False)
    hex_digest = hashlib.md5(df_norm_text.encode()).hexdigest()
    return hex_digest


def hash_primer_bed(bed_path: Path):
    """Hash a 7 column {scheme}.primer.bed file"""
    primer_bed_fields = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "poolName",
        "strand",
        "sequence",
    ]
    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=primer_bed_fields,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            strand=str,
            sequence=str,
        ),
    )
    df_norm = df.sort_values(["chromStart", "strand", "sequence"])
    concatenated_seqs = "|".join(df_norm["sequence"].tolist())
    chars_in_concatenated_seqs = set("".join(df_norm["sequence"].tolist()))
    allowed_chars = set("ACGT")
    if not chars_in_concatenated_seqs <= allowed_chars:
        raise RuntimeError("Illegal chars")
    hex_digest = hashlib.sha256(concatenated_seqs.encode()).hexdigest()
    return hex_digest


def hash_scheme_bed(bed_path: Path) -> str:
    """Hash a 6 column {scheme}.scheme.bed file"""
    pass


def hash_ref(ref_path: Path):
    record = SeqIO.read(ref_path, "fasta")
    hex_digest = hashlib.md5(str(record.seq).upper().encode()).hexdigest()
    return hex_digest


def build():
    """Build PHA4GE-compatible primer scheme bundle"""
    pass


def main():
    defopt.run(
        {
            "hash-bed": hash_bed,
            "hash-ref": hash_ref,
            "build": build,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()
