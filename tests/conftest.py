"""Shared synthetic input files; each test receives isolated temporary files."""
from pathlib import Path

import pytest


@pytest.fixture
def vcf(tmp_path: Path) -> Path:
    vcf = tmp_path / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC1-ExC1-xgenV1\tC2-ExC1-xgenV1\n"
        "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1/1\n"
    )
    return vcf


@pytest.fixture
def cohort(tmp_path: Path) -> Path:
    cohort = tmp_path / "cohort.csv"
    cohort.write_text("id,label\nC1-ExC1-xgenV1,one\nC1-ExC1-xgenV1,duplicate\n")
    return cohort


@pytest.fixture
def catalog(tmp_path: Path) -> Path:
    catalog = tmp_path / "catalog.csv"
    catalog.write_text("chr_pos,name,rsid\nchr1_1,TEST,rs123\n")
    return catalog
