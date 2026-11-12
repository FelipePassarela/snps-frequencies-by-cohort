import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

import usp_gene_correlation.cli as cli
from usp_gene_correlation.cli import app
from usp_gene_correlation.ordering import read_gene_order


@pytest.fixture
def catalog(tmp_path: Path) -> Path:
    catalog = tmp_path / "catalog.csv"
    catalog.write_text("chr_pos,name,rsid\nchr1_1,TEST,rs123\n")
    return catalog


@pytest.fixture
def gene_order(tmp_path: Path) -> Path:
    order = tmp_path / "gene_order.csv"
    order.write_text("name\nLERFS\nHK1\nABSENT\n")
    return order


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
    assert "rs123" in result.stdout


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
    assert "--gene-order" in result.output


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


@pytest.mark.parametrize("apply_order", [False, True])
@pytest.mark.parametrize("write_csv", [False, True])
def test_cli_gene_order_is_optional_for_stdout_and_csv(
    unordered_vcf: Path,
    cohort: Path,
    unordered_catalog: Path,
    gene_order: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    apply_order: bool,
    write_csv: bool,
) -> None:
    arguments = [
        str(unordered_vcf),
        str(cohort),
        "--annotations",
        str(unordered_catalog),
    ]
    if apply_order:
        arguments.extend(["--gene-order", str(gene_order)])
    else:

        def unexpected_order_read(path: Path) -> tuple[str, ...]:
            raise AssertionError(
                f"Gene order must not be read without --gene-order: {path}"
            )

        monkeypatch.setattr(cli, "read_gene_order", unexpected_order_read)
    output_dir = tmp_path / "ordered_results"
    if write_csv:
        arguments.extend(["-o", str(output_dir)])

    result = CliRunner().invoke(app, arguments)

    assert result.exit_code == 0, result.output
    expected = ["chr1_4", "chr1_3", "chr1_2", "chr1_5", "chr1_1", "chr1_3"]
    if apply_order:
        expected = ["chr1_1", "chr1_3", "chr1_2", "chr1_3", "chr1_4", "chr1_5"]
    if write_csv:
        actual = pd.read_csv(output_dir / "cohort.csv", keep_default_na=False)
        assert actual["variant"].tolist() == expected
        if apply_order:
            assert actual["name"].tolist() == ["LERFS", "HK1", "HK1", "HK1", "", ""]
            assert actual["rsid"].tolist() == ["rs1", "rs3", "rs2", "rs3", "", ""]
            assert actual["0/0"].tolist() == [1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    else:
        actual_variants = [
            line.split()[0]
            for line in result.stdout.splitlines()
            if line.startswith("chr1_")
        ]
        assert actual_variants == expected
        assert "ABSENT" not in result.stdout


def test_cli_reads_only_selected_order_once_for_multiple_cohorts(
    vcf: Path,
    cohort: Path,
    catalog: Path,
    gene_order: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Path] = []

    def record_order_read(path: Path) -> tuple[str, ...]:
        calls.append(path)
        return read_gene_order(path)

    monkeypatch.setattr(cli, "read_gene_order", record_order_read)
    monkeypatch.chdir(tmp_path)
    other = tmp_path / "other.csv"
    other.write_bytes(cohort.read_bytes())
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(
        app,
        [
            str(vcf),
            str(cohort),
            str(other),
            "--annotations",
            str(catalog),
            "--gene-order",
            str(gene_order),
            "-o",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls == [gene_order]
    assert sorted(path.name for path in output_dir.glob("*.csv")) == [
        "cohort.csv",
        "other.csv",
    ]


@pytest.mark.parametrize("invalid_order", ["missing", "directory", "duplicate"])
def test_cli_rejects_invalid_gene_order_before_writing_outputs(
    vcf: Path,
    cohort: Path,
    catalog: Path,
    tmp_path: Path,
    invalid_order: str,
) -> None:
    path = tmp_path / "bad_order.csv"
    if invalid_order == "directory":
        path.mkdir()
    elif invalid_order == "duplicate":
        path.write_text("name\nHK1\nHK1\n")
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(
        app,
        [
            str(vcf),
            str(cohort),
            "--annotations",
            str(catalog),
            "--gene-order",
            str(path),
            "-o",
            str(output_dir),
        ],
    )

    assert result.exit_code == 2, result.output
    assert "--gene-order" in result.output
    if invalid_order == "duplicate":
        assert "Duplicate gene" in result.output
    assert not output_dir.exists()
