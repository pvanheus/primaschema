import subprocess
from pathlib import Path

import pytest

data_dir = "tests/data"


def run(cmd, cwd="./"):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=data_dir, shell=True, check=True, text=True, capture_output=True
    )


def test_version():
    run("primaschema --version")


def test_hash_bed():
    run_cmd = run("primaschema hash-bed eden.scheme.bed")
    assert "4a9cc7082d23d32b" in run_cmd.stdout


def test_hash_ref():
    run_cmd = run("primaschema hash-ref eden.reference.fasta")
    assert "105c82802b67521950854a851fc6eefd" in run_cmd.stdout
