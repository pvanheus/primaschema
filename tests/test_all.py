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


def test_build():
    run_cmd = run("primaschema build sars-cov-2/eden/v1/eden.scheme.bed")
    assert "4a9cc7082d23d32b" in run_cmd.stdout
