import subprocess
import sys
from pathlib import Path
from typing import TypeAlias

import pandas as pd
import pytest
from typer.testing import CliRunner

from usp_gene_correlation.annotations import (
    CsvAnnotationSource,
    SnpAnnotation,
    annotate_frequencies,
)
from usp_gene_correlation.cli import app, main
from usp_gene_correlation.cohorts import get_cohort_ids
from usp_gene_correlation.genotypes import genotype_frequencies, genotypes_matrix

Inputs: TypeAlias = tuple[Path, Path, Path]


def test_annotations_preserve_values_order_and_unknowns(tmp_path: Path) -> None:
    original = pd.DataFrame(
        {"0/0": [0.25, 1.0, 0.5]}, index=["chr2_2", "1_1", "chr2_2"]
    )
    catalog = tmp_path / "annotations.csv"
    catalog.write_text("chr_pos,name,rsid\n2_2,GENE,rs1\n", encoding="utf-8")
    source = CsvAnnotationSource(catalog)
    with pytest.warns(UserWarning, match="1_1"):
        actual = annotate_frequencies(original, source)
    pd.testing.assert_frame_equal(actual.drop(columns=["name", "rsid"]), original)
    assert actual.name.tolist() == ["GENE", "", "GENE"]
    assert actual.rsid.tolist() == ["rs1", "", "rs1"]
    assert list(original.columns) == ["0/0"]


@pytest.mark.parametrize(
    "contents",
    [
        "position,name,rsid\nchr1_1,A,rs1\n",
        "chr_pos,name,rsid\nchr1_1,A,rs1\n1_1,B,rs2\n",
        "chr_pos,name,rsid\nchr1_0,A,rs1\n",
        "chr_pos,name,rsid\nchr1_1,A,bad\n",
        "chr_pos,name,rsid\nchr1_1,A\n",
    ],
)
def test_invalid_catalog(tmp_path: Path, contents: str) -> None:
    path = tmp_path / "snps_annotations.csv"
    path.write_text(contents)
    with pytest.raises(ValueError):
        CsvAnnotationSource(path)


def test_bundled_catalog(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    source = CsvAnnotationSource.from_default_catalog()
    assert source.lookup("9_62897604") == SnpAnnotation(
        "chr9_62897604", "LERFS", "rs200443822"
    )
    assert source.lookup("chr2_46818170").name == "LINC01118"
    assert len(source._records) == 44


def test_frequencies() -> None:
    frame = pd.DataFrame({"chr1_1": ["0/0", "1|0", "1/1", "0/.", "./."]})
    result = genotype_frequencies(frame)
    assert result.loc["chr1_1"].tolist() == pytest.approx([0.2, 0.2, 0.2, 0.4])
    assert genotype_frequencies(pd.DataFrame({"v": ["0/0"]})).iloc[0].tolist() == [
        1,
        0,
        0,
        0,
    ]


@pytest.mark.parametrize(
    "frame", [pd.DataFrame(), pd.DataFrame({"v": ["1/2"]}), pd.DataFrame({"v": ["1"]})]
)
def test_unsupported_or_empty_genotypes(frame: pd.DataFrame) -> None:
    with pytest.raises(ValueError):
        genotype_frequencies(frame)


@pytest.fixture
def inputs(tmp_path: Path) -> Inputs:
    vcf = tmp_path / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC1-ExC1-xgenV1\tC2-ExC1-xgenV1\n"
        "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1/1\n"
    )
    cohort = tmp_path / "cohort.csv"
    cohort.write_text("id,label\nC1-ExC1-xgenV1,one\nC1-ExC1-xgenV1,duplicate\n")
    catalog = tmp_path / "catalog.csv"
    catalog.write_text("chr_pos,name,rsid\nchr1_1,TEST,rs123\n")
    return vcf, cohort, catalog


def test_vcf_cohort_and_cli(inputs: Inputs, tmp_path: Path) -> None:
    vcf, cohort, catalog = inputs
    assert get_cohort_ids([cohort]).iloc[:, 0].dropna().tolist() == ["C1-ExC1-xgenV1"]
    assert genotypes_matrix(vcf, ["C1-ExC1-xgenV1"]).iloc[0, 0] == "0/1"
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


def test_injected_provider_and_stdout(
    inputs: Inputs, capsys: pytest.CaptureFixture[str]
) -> None:
    vcf, cohort, _ = inputs

    class ApiLikeSource:
        def lookup(self, position: str) -> SnpAnnotation:
            return SnpAnnotation(position, "INJECTED", "rs456")

    main([str(vcf), str(cohort)], annotation_source=ApiLikeSource())
    assert "INJECTED" in capsys.readouterr().out


def test_invalid_cohort(tmp_path: Path) -> None:
    path = tmp_path / "invalid.csv"
    path.write_text("name,value\nfoo,bar\n")
    with pytest.raises(ValueError, match="No matching"):
        get_cohort_ids([path])


@pytest.mark.parametrize("option", ["-o", "--output_dir", "--output-dir"])
def test_typer_output_aliases_and_multiple_cohorts(
    inputs: Inputs, tmp_path: Path, option: str
) -> None:
    vcf, cohort, catalog = inputs
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


def test_typer_invalid_input_path(inputs: Inputs, tmp_path: Path) -> None:
    _, cohort, _ = inputs
    out = tmp_path / "out"
    result = CliRunner().invoke(
        app, [str(tmp_path / "missing.vcf"), str(cohort), "-o", str(out)]
    )
    assert result.exit_code == 2
    assert not out.exists()
