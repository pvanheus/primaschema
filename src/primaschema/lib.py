import hashlib
import importlib.util
import io
import json
import logging
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

import jsonschema
import yaml

import altair as alt
import pandas as pd

from natsort import natsorted
from Bio import SeqIO


# from linkml.generators.pydanticgen import PydanticGenerator
from linkml.generators.pythongen import PythonGenerator
from linkml_runtime.utils.schemaview import SchemaView
from linkml.validators import JsonSchemaDataValidator

from primaschema import primer_scheme_schema_path, organisms_path, manifest_schema_path


SCHEME_BED_FIELDS = ["chrom", "chromStart", "chromEnd", "name", "poolName", "strand"]
PRIMER_BED_FIELDS = SCHEME_BED_FIELDS + ["sequence"]
POSITION_FIELDS = ["chromStart", "chromEnd"]

logger = logging.getLogger("primaschema")


def scan(path):
    """Recursively yield DirEntry objects"""
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            yield from scan(entry.path)
        else:
            yield entry


def import_class_from_path(file_path, class_name="PrimerScheme"):
    spec = importlib.util.spec_from_file_location(class_name, file_path)
    if spec is None:
        raise ImportError(f"Failed to load schema from {file_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, class_name)


# def get_primer_schemes_path():
#     """Locate primer-schemes repo root using environment variable"""
#     env_var = "PRIMER_SCHEMES_PATH"
#     if (
#         env_var not in os.environ
#         or not (
#             Path(os.environ[env_var]).resolve()
#             / Path("schema")
#             / Path("primer_scheme.yml")
#         ).exists()
#     ):
#         raise RuntimeError(
#             f'Invalid or unset environment variable {env_var} ({os.environ.get(env_var)}).\n\nSet {env_var} to the path of a local copy of the primer-schemes repo to proceed. For example, do `git clone https://github.com/pha4ge/primer-schemes` followed by `export {env_var}="/path/to/primer-schemes"`'
#         )
#     return Path(os.environ[env_var]).resolve()


def hash_string(string: str) -> str:
    """Normalise case, sorting, terminal spaces & return prefixed 64b of SHA256 hex"""
    checksum = hashlib.sha256(str(string).strip().upper().encode()).hexdigest()[:16]
    return f"primaschema:{checksum}"


def parse_scheme_bed(bed_path: Path) -> pd.DataFrame:
    """Parse a 6 column scheme.bed bed file"""
    return pd.read_csv(
        bed_path,
        sep="\t",
        names=SCHEME_BED_FIELDS,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=int,
            strand=str,
        ),
    )


def parse_primer_bed(bed_path: Path) -> pd.DataFrame:
    """Parse a 7 column primer.bed bed file"""
    return pd.read_csv(
        bed_path,
        sep="\t",
        names=PRIMER_BED_FIELDS,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=int,
            strand=str,
            sequence=str,
        ),
    )


def normalise_primer_bed_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    - Removes terminal whitespace
    - Normalises case
    - Sorts by chromStart, chromEnd, poolName, strand, sequence
    - Removes duplicate records, collapsing alts with same coords if backfilled from ref
    """
    df["sequence"] = df["sequence"].str.strip().str.upper()
    df = df.sort_values(
        ["chromStart", "chromEnd", "poolName", "strand", "sequence"]
    ).drop_duplicates()
    return df


def hash_primer_bed_df(df: pd.DataFrame) -> str:
    """
    Returns prefixed SHA256 digest from stringified dataframe
    """
    string = df[["chromStart", "chromEnd", "poolName", "strand", "sequence"]].to_csv(
        index=False
    )
    return hash_string(string)


def hash_primer_bed(bed_path: Path):
    """Hash a 7 column primer.bed file"""
    df = parse_primer_bed(bed_path)
    return hash_primer_bed_df(df)


def hash_scheme_bed(bed_path: Path, fasta_path: Path) -> str:
    """
    Hash a 6 column scheme.bed file by first converting to 7 column primer.bed
    """
    logger.info("Hashing scheme.bed using reference backfill")
    ref_record = SeqIO.read(fasta_path, "fasta")
    df = parse_scheme_bed(bed_path)
    records = df.to_dict("records")
    for r in records:
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            r["sequence"] = str(ref_record.seq[start_pos:end_pos])
        elif r["strand"] == "-":
            r["sequence"] = str(ref_record.seq[start_pos:end_pos].reverse_complement())
        else:
            raise RuntimeError(f"Invalid strand for BED record {r}")
    bed7_df = pd.DataFrame(records)
    return hash_primer_bed_df(bed7_df)


def convert_primer_bed_to_scheme_bed(bed_path: Path):
    df = parse_primer_bed(bed_path).drop("sequence", axis=1)
    with io.StringIO() as csv_buffer:
        df.to_csv(csv_buffer, sep="\t", header=False, index=False)
        bed_str = csv_buffer.getvalue()
    return bed_str


def convert_scheme_bed_to_primer_bed(bed_path: Path, fasta_path: Path) -> str:
    ids_seqs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    df = parse_scheme_bed(bed_path)
    records = df.to_dict("records")
    for r in records:
        chrom = r["chrom"].partition(" ")[0]  # Use chrom name before first space
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            r["sequence"] = str(ids_seqs[chrom].seq[start_pos:end_pos])
        else:
            r["sequence"] = str(
                ids_seqs[chrom].seq[start_pos:end_pos].reverse_complement()
            )
    df = pd.DataFrame(records)
    with io.StringIO() as csv_buffer:
        df.to_csv(csv_buffer, sep="\t", header=False, index=False)
        bed_str = csv_buffer.getvalue()
    return bed_str


def hash_bed(bed_path: Path) -> str:
    bed_type = infer_bed_type(bed_path)
    if bed_type == "primer":
        checksum = hash_primer_bed(bed_path)
    else:  # bed_type == "scheme"
        checksum = hash_scheme_bed(
            bed_path=bed_path, fasta_path=bed_path.parent / "reference.fasta"
        )
    return checksum


def hash_ref(ref_path: Path):
    record = SeqIO.read(ref_path, "fasta")
    return hash_string(record.seq)


def count_tsv_columns(bed_path: Path) -> int:
    return len(pd.read_csv(bed_path, sep="\t").columns)


def parse_yaml(path) -> dict:
    with open(path, "r") as fh:
        return yaml.safe_load(fh)


def validate_yaml_with_json_schema(yaml_path: Path, schema_path: Path):
    yaml_data = parse_yaml(yaml_path)
    with open(schema_path, "r") as schema_fh:
        schema = json.load(schema_fh)
    return jsonschema.validate(yaml_data, schema=schema)


def validate_with_linkml_schema(yaml_path: Path, schema_path: Path, full: bool = False):
    schema_view = SchemaView(schema_path)
    schema_gen = PythonGenerator(schema_view.schema)
    schema_compiled = schema_gen.compile_module()
    data = parse_yaml(yaml_path)
    data_instance = schema_compiled.PrimerScheme(**data)
    # print(yaml_dumper.dumps(data_instance))
    validator = JsonSchemaDataValidator(schema_view.schema)
    validator.validate_object(data_instance)


# def validate_with_linkml_schema(yaml_path: Path, full: bool = False):
#     data = parse_yaml(yaml_path)
#     schema_path = get_primer_schemes_path() / "schema/primer_scheme.yml"
#     pythonised_schema_path = get_primer_schemes_path() / "schema/primer_scheme.py"
#     if full:
#         schema_view = SchemaView(schema_path)
#         schema_gen = PythonGenerator(schema_view.schema)
#         schema_compiled = schema_gen.compile_module()
#         schema_compiled.PrimerScheme(**data)  # Errors on validation failure
#     else:
#         if not pythonised_schema_path.exists():
#             run(f"gen-python {schema_path} > {pythonised_schema_path}")
#             logger.info(f"Wrote Pythonised schema to {pythonised_schema_path}")
#             print(run("ls").stdout)
#         PrimerScheme = import_class_from_path(pythonised_schema_path)
#         PrimerScheme(**data)  # Errors on validation failure


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
        logger.info(f"Detected {bed_type} bed file with {bed_columns} columns")

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


def validate(scheme_dir: Path, full: bool = False, force: bool = False):
    logger.info(f"Validating {scheme_dir}")
    validate_bed(scheme_dir / "primer.bed", bed_type="primer")
    validate_with_linkml_schema(
        yaml_path=scheme_dir / "info.yml", schema_path=primer_scheme_schema_path
    )
    # validate_with_linkml_schema(yaml_path=scheme_dir / "info.yml", full=full)
    scheme = parse_yaml(scheme_dir / "info.yml")
    existing_primer_checksum = scheme.get("primer_checksum")
    existing_reference_checksum = scheme.get("reference_checksum")
    primer_checksum = hash_bed(scheme_dir / "primer.bed")
    reference_checksum = hash_ref(scheme_dir / "reference.fasta")
    if (
        existing_primer_checksum
        and not primer_checksum == existing_primer_checksum
        and not force
    ):
        raise RuntimeError(
            f"Calculated and documented primer checksums do not match ({primer_checksum} and {existing_primer_checksum})"
        )
    elif not primer_checksum == existing_primer_checksum:
        logging.warning(
            f"Calculated and documented primer checksums do not match ({primer_checksum} and {existing_primer_checksum})"
        )
    if (
        existing_reference_checksum
        and not reference_checksum == existing_reference_checksum
        and not force
    ):
        raise RuntimeError(
            f"Calculated and documented reference checksums do not match ({reference_checksum} and {existing_reference_checksum})"
        )
    elif not reference_checksum == existing_reference_checksum:
        logging.warning(
            f"Calculated and documented reference checksums do not match ({reference_checksum} and {existing_reference_checksum})"
        )
    logger.info(f"Validation successful for {scheme.get('name')} ")


def validate_recursive(root_dir: Path, full: bool = False, force: bool = False):
    """Validate all schemes in a directory tree"""
    schemes_paths = {}
    for entry in scan(root_dir):
        if entry.is_file() and entry.name == "info.yml":
            scheme_info = parse_yaml(entry.path)
            scheme_name = scheme_info["name"]
            scheme_dir = Path(entry.path).parent
            schemes_paths[scheme_name] = scheme_dir

    for scheme_name, path in schemes_paths.items():
        validate(scheme_dir=path, full=full, force=force)


def build(
    scheme_dir: Path,
    out_dir: Path = Path(),
    full: bool = False,
    force: bool = False,
    nested: bool = True,
):
    """
    Build a PHA4GE primer scheme bundle.
    Given a directory path containing info.yml, reference.fasta, and either
    primer.bed or reference.bed, generate a directory containing info.yml including
    primer and reference checksums and a canonical primer.bed representation.
    """
    validate(scheme_dir=scheme_dir, full=full, force=force)
    scheme = parse_yaml(scheme_dir / "info.yml")
    if nested:
        family = Path(scheme["name"].partition("-")[0])
        version = Path(scheme["name"].partition("-")[2])
        out_dir = Path("built") / scheme["organism"] / family / version
    else:
        out_dir = Path("built") / scheme["name"]
    try:
        out_dir.mkdir(parents=True, exist_ok=force)
    except FileExistsError:
        raise FileExistsError(f"Output directory {out_dir} already exists")
    scheme["primer_checksum"] = hash_bed(scheme_dir / "primer.bed")
    scheme["reference_checksum"] = hash_ref(scheme_dir / "reference.fasta")
    with open(out_dir / "info.yml", "w") as scheme_fh:
        logger.info(f"Writing info.yml to {out_dir}/info.yml")
        yaml.dump(scheme, scheme_fh, sort_keys=False)
    logger.info(f"Copying primer.bed to {out_dir}/primer.bed")
    shutil.copy(scheme_dir / "primer.bed", out_dir)
    logger.info(f"Copying reference.fasta to {out_dir}/reference.fasta")
    shutil.copy(scheme_dir / "reference.fasta", out_dir)
    logger.info(f"Writing scheme.bed to {out_dir}/scheme.bed")
    scheme_bed_str = convert_primer_bed_to_scheme_bed(bed_path=out_dir / "primer.bed")
    with open("scheme.bed", "w") as fh:
        fh.write(scheme_bed_str)
    shutil.copy("scheme.bed", out_dir.resolve())
    os.remove("scheme.bed")


def build_recursive(
    root_dir: Path, full: bool = False, force: bool = False, nested: bool = False
):
    """Build all schemes in a directory tree"""
    schemes_paths = {}
    for entry in scan(root_dir):
        if entry.is_file() and entry.name == "info.yml":
            scheme = parse_yaml(entry.path)
            scheme_dir = Path(entry.path).parent
            schemes_paths[scheme.get("name")] = scheme_dir
    for scheme, path in schemes_paths.items():
        build(scheme_dir=path, full=full, force=force)


def build_manifest(root_dir: Path, out_dir: Path = Path()):
    """Build manifest of schemes inside the specified directory"""
    organisms = parse_yaml(organisms_path)
    manifest = {
        "schema_version": "0.9.0",
        "metadata": "The PHA4GE list of tiling amplicon primer schemes",
        "repository": "https://github.com/pha4ge/primer-schemes",
        "latest_doi": "",
        "license": "CC-BY-4.0",
        "organisms": organisms,
    }

    names_schemes = {}
    families_names = defaultdict(list)
    for entry in scan(root_dir):
        if entry.is_file() and entry.name == "info.yml":
            scheme = parse_yaml(entry.path)
            name = scheme["name"]
            names_schemes[name] = scheme
            family, _, version = scheme["name"].partition("-")
            families_names[family].append(name)

    families_data = []
    for family, names in sorted(families_names.items()):
        family_data = {}
        family_data["family"] = family
        family_example_name = families_names[family][0]
        family_data["organism"] = names_schemes[family_example_name]["organism"]
        versions_data = []
        definition_url = names_schemes[name].get("definition_url", "")
        for name in sorted(names):
            if names_schemes[name].get("display_name"):
                display_name = names_schemes[name]["display_name"]
            else:
                display_name = name
            versions_data.append(
                {
                    "name": name,
                    "display_name": display_name,
                    "version": name.partition("-")[2],
                    "repository": definition_url,
                }
            )
            logger.info(f"Reading {name}")
        family_data["versions"] = versions_data
        families_data.append(family_data)
    manifest["schemes"] = families_data

    manifest_file_name = "index.yml"
    with open(out_dir / manifest_file_name, "w") as fh:
        logger.info(f"Writing {manifest_file_name} to {out_dir}/{manifest_file_name}")
        yaml.dump(data=manifest, stream=fh, sort_keys=False)
    validate_yaml_with_json_schema(
        yaml_path=out_dir / manifest_file_name, schema_path=manifest_schema_path
    )


def diff(bed1_path: Path, bed2_path: Path, only_positions: bool = False):
    """Show symmetric differences between records in two primer.bed files"""
    df1 = parse_primer_bed(bed1_path).assign(origin="bed1")
    df2 = parse_primer_bed(bed2_path).assign(origin="bed2")
    if only_positions:
        column_subset = POSITION_FIELDS
    else:
        column_subset = PRIMER_BED_FIELDS
    return pd.concat([df1, df2]).drop_duplicates(subset=column_subset, keep=False)


def show_non_ref_alts(scheme_dir: Path):
    """Show primer records with sequences not matching the reference sequence"""
    bed_path = scheme_dir / "primer.bed"
    with TemporaryDirectory() as temp_dir:
        convert_scheme_bed_to_primer_bed(
            bed_path=scheme_dir / "scheme.bed",
            fasta_path=scheme_dir / "reference.fasta",
        )
        return diff(bed1_path=bed_path, bed2_path=Path(temp_dir) / "primer.bed")


def compute_intervals(bed_path: Path) -> dict[str, dict[str, (int, int)]]:
    """
    find primer positions for all primers in the bed file and compute maximum
    interval between primers of the same name
    """
    primer_name_re = re.compile(r"^(?P<name>.*)_(LEFT|RIGHT)(_.+)?$")
    eden_primer_name_re = re.compile(r"^(?P<name>.*_[AB][0-9])(F|R)_\d+$")
    all_intervals: dict[str, dict[str, (int, int)]] = {}
    for line in open(bed_path):
        line_parts = line.strip().split("\t")
        if len(line_parts) < 6:
            # skip lines that don't have at least 6 fields
            continue
        chrom, start, end, name, _, strand = line.strip().split("\t")[:6]
        if chrom not in all_intervals:
            all_intervals[chrom] = {}
        intervals = all_intervals[chrom]
        primer_match = primer_name_re.match(name)
        if not primer_match:
            # the Eden scheme has a unique primer name format
            primer_match = eden_primer_name_re.match(name)
            if not primer_name_re:
                raise ValueError(f"Invalid primer name {name}")
        primer_name = primer_match.group("name")
        if strand == "+":
            start_pos = int(start)
            end_pos = -1
        if strand == "-":
            start_pos = sys.maxsize
            end_pos = int(end)
        prev_start, prev_end = intervals.get(primer_name, (sys.maxsize, -1))
        intervals[primer_name] = (min(prev_start, start_pos), max(prev_end, end_pos))
    return all_intervals


def plot(
    bed_path: Path, amplicon_name_index: int = 1, out_path: Path = Path("plot.html")
) -> None:
    """
    Plot amplicon and primer positions from a 7 column primer.bed file
    Requires primers to be named {scheme-name}_{amplicon-number}…
    Plots one row per reference chromosome
    Supported out_path extensions: html (interactive), pdf, png, svg
    """
    bed_df = parse_primer_bed(bed_path)
    bed_df["amplicon"] = bed_df["name"].str.split("_").str[amplicon_name_index]
    logging.debug(pd.concat([bed_df.head(4), bed_df.tail(4)]))
    amp_df = (
        bed_df.groupby(["chrom", "amplicon"])
        .agg(min_start=("chromStart", "min"), max_end=("chromEnd", "max"))
        .reset_index()
    )
    amp_df["is_amplicon"] = True
    sorted_amplicons = natsorted(bed_df["amplicon"].unique())
    logging.debug(amp_df)

    bed_df["is_amplicon"] = False
    amp_df = amp_df.rename(columns={"min_start": "chromStart", "max_end": "chromEnd"})

    combined_df = pd.concat([bed_df, amp_df], ignore_index=True)

    primer_marks = (
        alt.Chart(combined_df)
        .transform_filter(alt.datum.is_amplicon == False)  # noqa
        .mark_line(size=15)
        .encode(
            x=alt.X("chromStart:Q", title=None),
            x2="chromEnd:Q",
            y=alt.Y("amplicon:O", sort=sorted_amplicons, scale=alt.Scale(padding=0)),
            color=alt.Color("strand:N").scale(scheme="set2"),
            tooltip=[
                alt.Tooltip("name:N", title="Primer name"),
                alt.Tooltip("chromStart:Q", title="start"),
                alt.Tooltip("chromEnd:Q", title="end"),
            ],
        )
        .properties(
            width=1000,
        )
    )
    amplicon_marks = (
        alt.Chart(combined_df)
        .transform_filter(alt.datum.is_amplicon == True)  # noqa
        .mark_rule(size=2)
        .encode(
            x=alt.X("chromStart:Q", title=None),
            x2="chromEnd:Q",
            y=alt.Y("amplicon:O", sort=sorted_amplicons, scale=alt.Scale(padding=0)),
            tooltip=[
                alt.Tooltip("amplicon:N", title="Amplicon name"),
                alt.Tooltip("chromStart:Q", title="Min primer start"),
                alt.Tooltip("chromEnd:Q", title="Max primer end"),
            ],
        )
        .properties(
            width=1000,
        )
    )
    combined_chart = alt.layer(primer_marks, amplicon_marks).facet(
        row=alt.Row("chrom:O", header=alt.Header(labelOrient="top"), title="")
    )
    combined_chart.interactive().save(str(out_path))
    logger.info(f"Plot saved ({out_path})")
