import hashlib
import json
import logging
import shutil
from pathlib import Path
from typing import Literal

import jsonschema
import pandas as pd
import yaml
from Bio import SeqIO

from primaschema import schema_dir

SCHEME_BED_FIELDS = ["chrom", "chromStart", "chromEnd", "name", "poolName", "strand"]
PRIMER_BED_FIELDS = SCHEME_BED_FIELDS + ["sequence"]


def hash_sequences(sequences: list[str]) -> str:
    """
    Normalise by uppercasing and sorting sequences
    Check for non-ACGT chars
    Return the sha256 hex digest of the comma delimited normalised sequence list
    """
    sequences = sorted(map(str.upper, sequences))
    concatenated_seqs = ",".join(sequences)
    chars_in_concatenated_seqs = set("".join(sequences))
    allowed_chars = set("ACGT")
    if not chars_in_concatenated_seqs <= allowed_chars:
        raise RuntimeError("Found illegal characters in primer sequences")
    hex_digest = hashlib.sha256(concatenated_seqs.encode()).hexdigest()
    return f"primaschema:{hex_digest}"


def hash_primer_bed(bed_path: Path):
    """Hash a 7 column {scheme}.primer.bed file"""
    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=PRIMER_BED_FIELDS,
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
    return hash_sequences(df["sequence"].tolist())


def hash_scheme_bed(bed_path: Path, fasta_path: Path) -> str:
    """
    Hash a 6 column {scheme}.scheme.bed file

    Seqs need to first be retrieved and sorted against in order to normalise prior to hashing
    """
    ref_record = SeqIO.read(fasta_path, "fasta")
    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=SCHEME_BED_FIELDS,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            strand=str,
        ),
    )
    primer_sequences = []
    for r in df.to_dict("records"):
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            primer_sequences.append(str(ref_record.seq[start_pos:end_pos]))
        else:
            primer_sequences.append(
                str(ref_record.seq[start_pos:end_pos].reverse_complement())
            )
    return hash_sequences(primer_sequences)


def convert_scheme_bed_to_primer_bed(
    bed_path: Path, fasta_path: Path, out_dir: Path = Path(), force: bool = False
):
    ref_record = SeqIO.read(fasta_path, "fasta")
    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=SCHEME_BED_FIELDS,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            strand=str,
        ),
    )
    records = df.to_dict("records")
    for r in records:
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            r["sequence"] = str(ref_record.seq[start_pos:end_pos])
        else:
            r["sequence"] = str(ref_record.seq[start_pos:end_pos].reverse_complement())
    df = pd.DataFrame(records)
    df.to_csv(out_dir / "primer.bed", sep="\t", header=False, index=False)


def hash_bed(bed_path: Path) -> str:
    bed_type = infer_bed_type(bed_path)
    if bed_type == "primer":
        checksum = hash_primer_bed(bed_path)
    else:  # primary_bed_type == "scheme"
        checksum = hash_scheme_bed(
            bed_path=bed_path, fasta_path=bed_path.parent / "reference.fasta"
        )
    return checksum


def hash_ref(ref_path: Path):
    record = SeqIO.read(ref_path, "fasta")
    checksum = hashlib.sha256(str(record.seq).upper().encode()).hexdigest()
    return f"primaschema:{checksum}"


def count_tsv_columns(bed_path: Path) -> int:
    return len(pd.read_csv(bed_path, sep="\t").columns)


def parse_scheme(scheme_path) -> dict:
    with open(scheme_path, "r") as scheme_fh:
        return yaml.safe_load(scheme_fh)


def validate_yaml(scheme_path):
    schema_path = schema_dir / "scheme_schema.json"
    with open(schema_path, "r") as schema_fh:
        schema = json.load(schema_fh)
    scheme = parse_scheme(scheme_path)
    return jsonschema.validate(scheme, schema=schema)


def validate_bed(bed_path: Path, bed_type=Literal["primer", "scheme"]):
    bed_columns = count_tsv_columns(bed_path)
    if bed_type == "primer" and bed_columns != 7:
        raise RuntimeError(
            f"Primer bed files should have 7 columns: {PRIMER_BED_FIELDS}"
        )
    elif bed_type == "scheme" and bed_columns != 6:
        raise RuntimeError(
            f"Scheme bed files should have 6 columns: {SCHEME_BED_FIELDS}"
        )
    else:
        logging.info(
            f"{bed_type.title()} bed file has the expected number of columns ({bed_columns})"
        )

    if bed_type == "primer":
        hash_primer_bed(bed_path)
    elif bed_type == "scheme":
        hash_scheme_bed(
            bed_path=bed_path, fasta_path=bed_path.parent / "reference.fasta"
        )


def infer_bed_type(bed_path: Path) -> str:
    bed_columns = count_tsv_columns(bed_path)
    if bed_columns == 7:
        bed_type = "primer"
    elif bed_columns == 6:
        bed_type = "scheme"
    else:
        raise RuntimeError(
            "Bed file shoud have either 6 columns (scheme.bed) or 7 column (primer.bed)"
        )
    return bed_type


def infer_primary_bed_type(scheme_dir: Path) -> str:
    if (scheme_dir / "primer.bed").exists():
        primary_bed_type = "primer"
    elif (scheme_dir / "scheme.bed").exists():
        primary_bed_type = "scheme"
    else:
        raise RuntimeError("Failed to discover a bed file")
    return primary_bed_type


def validate(scheme_dir: Path):
    primary_bed_type = infer_primary_bed_type(scheme_dir)
    validate_bed(scheme_dir / f"{primary_bed_type}.bed", bed_type=primary_bed_type)
    validate_yaml(scheme_dir / "info.yaml")
    scheme = parse_scheme(scheme_dir / "info.yaml")
    existing_primer_checksum = scheme.get("primer_checksum")
    existing_reference_checksum = scheme.get("reference_checksum")
    primer_checksum = hash_bed(scheme_dir / f"{primary_bed_type}.bed")
    reference_checksum = hash_ref(scheme_dir / "reference.fasta")
    if existing_primer_checksum and not primer_checksum == existing_primer_checksum:
        raise RuntimeError(
            f"Calculated and documented primer checksums do not match ({primer_checksum} and {existing_primer_checksum})"
        )
    if (
        existing_reference_checksum
        and not reference_checksum == existing_reference_checksum
    ):
        raise RuntimeError(
            f"Calculated and documented reference checksums do not match ({reference_checksum} and {existing_reference_checksum})"
        )


def build(scheme_dir: Path, out_dir: Path = Path(), force: bool = False):
    """
    Build a PHA4GE primer scheme bundle.
    Given a directory path containing info.yaml, reference.fasta, and either
    primer.bed or reference.bed, generate a directory containing info.yaml including
    primer and reference checksums and a canonical primer.bed representation.
    """
    validate(scheme_dir)
    primary_bed_type = infer_primary_bed_type(scheme_dir)
    scheme = parse_scheme(scheme_dir / "info.yaml")
    from pprint import pprint

    out_dir = Path(scheme["name"])
    try:
        out_dir.mkdir(exist_ok=force)
    except FileExistsError:
        raise FileExistsError(f"Output directory {out_dir} already exists")
    if not scheme.get("primer_checksum"):
        scheme["primer_checksum"] = hash_bed(scheme_dir / f"{primary_bed_type}.bed")
    if not scheme.get("reference_checksum"):
        scheme["reference_checksum"] = hash_ref(scheme_dir / "reference.fasta")
    with open(out_dir / "info.yaml", "w") as scheme_fh:
        logging.info(f"Writing {out_dir}/info.yaml")
        yaml.dump(scheme, scheme_fh, sort_keys=False)
    if primary_bed_type == "scheme":
        logging.info("Generating primer.bed from scheme.bed and reference.fasta")
        convert_scheme_bed_to_primer_bed(
            bed_path=scheme_dir / "scheme.bed",
            fasta_path=scheme_dir / "reference.fasta",
            out_dir=out_dir,
            force=force,
        )
    else:
        logging.info(f"Copying primer.bed to {out_dir}/primer.bed")
        shutil.copy(scheme_dir / "primer.bed", out_dir)
    logging.info(f"Copying reference.fasta to {out_dir}/reference.fasta")
    shutil.copy(scheme_dir / "reference.fasta", out_dir)
