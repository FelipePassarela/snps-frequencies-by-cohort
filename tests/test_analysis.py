from pathlib import Path

import pandas as pd
import pytest

from usp_gene_correlation.analysis import calculate_cohort_frequencies
from usp_gene_correlation.annotations import (
    CsvAnnotationSource,
    SnpAnnotation,
    SnpAnnotationSource,
)
from usp_gene_correlation.cohort import Cohort


class ApiLikeSource(SnpAnnotationSource):
    def lookup(self, chr_pos: str) -> SnpAnnotation:
        return SnpAnnotation(chr_pos, "INJECTED", "rs456")


def test_processing_accepts_injected_source(vcf: Path) -> None:
    result = calculate_cohort_frequencies(
        vcf, Cohort("test", ("C1-ExC1-xgenV1",)), ApiLikeSource()
    )
    assert result.loc["chr1_1", "name"] == "INJECTED"
    assert result.loc["chr1_1", "rsid"] == "rs456"
    assert result.loc["chr1_1", "0/1"] == 1


@pytest.mark.parametrize("gene_order", [None, ("LERFS", "HK1", "ABSENT")])
def test_processing_preserves_values_and_only_sorts_when_requested(
    unordered_vcf: Path,
    unordered_catalog: Path,
    gene_order: tuple[str, ...] | None,
) -> None:
    result = calculate_cohort_frequencies(
        unordered_vcf,
        Cohort("test", ("C1-ExC1-xgenV1",)),
        CsvAnnotationSource(unordered_catalog),
        gene_order=gene_order,
    )
    expected = pd.DataFrame(
        [
            ["", "", 0.0, 1.0, 0.0, 0.0],
            ["HK1", "rs3", 0.0, 0.0, 1.0, 0.0],
            ["HK1", "rs2", 0.0, 1.0, 0.0, 0.0],
            ["", "", 0.0, 0.0, 0.0, 1.0],
            ["LERFS", "rs1", 1.0, 0.0, 0.0, 0.0],
            ["HK1", "rs3", 1.0, 0.0, 0.0, 0.0],
        ],
        index=pd.Index(
            ["chr1_4", "chr1_3", "chr1_2", "chr1_5", "chr1_1", "chr1_3"],
            name="variant",
        ),
        columns=pd.Index(["name", "rsid", "0/0", "0/1", "1/1", "./."], name="genotype"),
    )
    if gene_order is not None:
        expected = expected.iloc[[4, 1, 2, 5, 0, 3]]

    pd.testing.assert_frame_equal(result, expected)
