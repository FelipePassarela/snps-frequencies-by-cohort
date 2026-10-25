import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from usp_gene_correlation.cli import app


def test_module_cli(vcf: Path, cohort: Path, catalog: Path, tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "usp_gene_correlation",
            str(vcf),
            str(cohort),
            "--annotations",
            str(catalog),
            "-o",
            str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    df = pd.read_csv(out / "cohort.csv")
    assert df.columns.tolist() == [
        "variant",
        "name",
        "rsid",
        "0/0",
        "0/1",
        "1/1",
        "./.",
    ]
    assert df.loc[0, "name"] == "TEST"
    assert df.loc[0, "0/1"] == 1


def test_stdout(vcf: Path, cohort: Path, catalog: Path) -> None:
    result = CliRunner().invoke(
        app, [str(vcf), str(cohort), "--annotations", str(catalog)]
    )
    assert result.exit_code == 0, result.output
    assert "TEST" in result.output
    assert "rs123" in result.output


@pytest.mark.parametrize("option", ["-o", "--output_dir", "--output-dir"])
def test_typer_output_aliases_and_multiple_cohorts(
    vcf: Path, cohort: Path, catalog: Path, tmp_path: Path, option: str
) -> None:
    other = tmp_path / "other.csv"
    other.write_bytes(cohort.read_bytes())
    out = tmp_path / "results"
    result = CliRunner().invoke(
        app,
        [
            str(vcf),
            str(cohort),
            str(other),
            "--annotations",
            str(catalog),
            option,
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert sorted(p.name for p in out.glob("*.csv")) == ["cohort.csv", "other.csv"]
    assert pd.read_csv(out / "other.csv").loc[0, "rsid"] == "rs123"


def test_typer_help() -> None:
    result = CliRunner().invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "--annotations" in result.output


def test_typer_missing_arguments() -> None:
    result = CliRunner().invoke(app, [])
    assert result.exit_code == 2


def test_typer_invalid_input_path(cohort: Path, tmp_path: Path) -> None:
    out = tmp_path / "out"
    result = CliRunner().invoke(
        app, [str(tmp_path / "missing.vcf"), str(cohort), "-o", str(out)]
    )
    assert result.exit_code == 2
    assert not out.exists()
