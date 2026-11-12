"""Shared synthetic inputs; each test receives isolated files and data."""

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
def unordered_vcf(vcf: Path) -> Path:
    header = "\n".join(vcf.read_text().splitlines()[:-1]) + "\n"
    vcf.write_text(
        header + "chr1\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t./.\n"
        "chr1\t3\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t1/1\n"
        "chr1\t2\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n"
        "chr1\t5\t.\tA\tG\t.\tPASS\t.\tGT\t./.\t0/0\n"
        "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"
        "chr1\t3\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\n"
    )
    return vcf


@pytest.fixture
def unordered_catalog(tmp_path: Path) -> Path:
    catalog = tmp_path / "unordered_catalog.csv"
    catalog.write_text(
        "chr_pos,name,rsid\n"
        "1_2,HK1,rs2\n"
        "chr1_6,ABSENT,rs6\n"
        "chr1_1,LERFS,rs1\n"
        "chr1_3,HK1,rs3\n"
    )
    return catalog
