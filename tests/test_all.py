import subprocess
from pathlib import Path

import pytest

import primaschema.lib as lib

data_dir = Path("tests/data")


def run(cmd, cwd=data_dir):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def test_cli_version():
    run("primaschema --version")


def test_hash_ref():
    assert (
        lib.hash_ref("tests/data/primer-schemes/eden/v1/reference.fasta")
        == "primaschema:7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"
    )


def test_cli_hash_ref():
    run_cmd = run("primaschema hash-ref primer-schemes/eden/v1/reference.fasta")
    assert (
        "primaschema:7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"
        in run_cmd.stdout
    )


def test_cli_hash_primer_bed():
    run_cmd = run("primaschema hash-bed primer-schemes/artic/v4.1/primer.bed")
    assert (
        "primaschema:bca2e6513bc5a45d664f7ddb43d6b63f3a67509b4381641ee030d333d29569c5"
        in run_cmd.stdout
    )


def test_cli_scheme_bed():
    run_cmd = run("primaschema hash-bed primer-schemes/artic/v4.1/scheme.bed")
    assert (
        "primaschema:bca2e6513bc5a45d664f7ddb43d6b63f3a67509b4381641ee030d333d29569c5"
        in run_cmd.stdout
    )


def test_artic_v41_scheme_hash_matches_primer_hash():
    scheme_bed_hash = lib.hash_scheme_bed(
        "tests/data/primer-schemes/artic/v4.1/scheme.bed",
        "tests/data/primer-schemes/artic/v4.1/reference.fasta",
    )
    primer_bed_hash = lib.hash_primer_bed(
        "tests/data/primer-schemes/artic/v4.1/primer.bed"
    )
    assert scheme_bed_hash == primer_bed_hash


def test_eden_v1_schema():
    lib.validate_yaml(f"{data_dir}/primer-schemes/eden/v1/info.yaml")


def test_artic_v41_schema():
    lib.validate_yaml(f"{data_dir}/primer-schemes/artic/v4.1/info.yaml")


# Needs updating since reverting hash function to consume coordinates again. Needs BEDs creating for this case
# def test_checksum_case_normalisation():
#     seqs_a = ["ACGT", "CAGT"]
#     seqs_b = ["ACGT", "cagt"]
#     assert lib.hash_sequences(seqs_a) == lib.hash_sequences(seqs_b)
def test_checksum_case_normalisation():
    assert lib.hash_bed(
        data_dir / "broken/different-case/eden-v1.primer.bed"
    ) == lib.hash_bed(data_dir / "broken/different-case/eden-v1-modified.primer.bed")


def test_validate_artic_v41():
    lib.validate(data_dir / "primer-schemes/artic/v4.1")


def test_validate_fail_five_columns():
    with pytest.raises(RuntimeError):
        lib.validate(data_dir / "broken/five-columns")


def test_hash_bed():
    lib.hash_bed(data_dir / "primer-schemes/artic/v4.1/primer.bed")
    lib.hash_bed(data_dir / "primer-schemes/artic/v4.1/scheme.bed")


def test_build_from_primer_bed():
    run_cmd = run("primaschema build primer-schemes/artic/v4.1 --force")
    run("rm -rf artic-v4.1")


def test_build_from_scheme_bed():
    run_cmd = run("primaschema build primer-schemes/eden/v1 --force")
    run("rm -rf eden-v1")


def test_build_recursively():
    lib.build_recursively(data_dir / "primer-schemes", force=True)
    run("rm -rf built", cwd="./")


def test_diff():
    run_cmd = run(
        "primaschema diff primer-schemes/midnight/v1/primer.bed primer-schemes/midnight/v2/primer.bed"
    )
    assert (
        """chrom  chromStart  chromEnd                      name poolName strand                 sequence origin
MN908947.3       27784     27808 SARS-CoV-2_28_LEFT_27837T        2      + TTTGTGCTTTTTAGCCTTTCTGTT   bed2"""
        in run_cmd.stdout
    )
