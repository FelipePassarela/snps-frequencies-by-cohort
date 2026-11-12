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
from usp_gene_correlation.ordering import FrequencyOrdering, GeneOrder, PreserveOrder


class ApiLikeSource(SnpAnnotationSource):
    def lookup(self, chr_pos: str) -> SnpAnnotation:
        return SnpAnnotation(chr_pos, "INJECTED", "rs456")


class ReverseOrder:
    def apply(self, frequencies: pd.DataFrame) -> pd.DataFrame:
        return frequencies.iloc[::-1]


def test_processing_accepts_injected_source(vcf: Path) -> None:
    result = calculate_cohort_frequencies(
        vcf, Cohort("test", ("C1-ExC1-xgenV1",)), ApiLikeSource()
    )
    assert result.loc["chr1_1", "name"] == "INJECTED"
    assert result.loc["chr1_1", "rsid"] == "rs456"
    assert result.loc["chr1_1", "0/1"] == 1


@pytest.mark.parametrize(
    ("ordering", "expected_rows"),
    [
        (PreserveOrder(), [0, 1, 2, 3, 4, 5]),
        (GeneOrder(("LERFS", "HK1", "ABSENT")), [4, 1, 2, 5, 0, 3]),
        (ReverseOrder(), [5, 4, 3, 2, 1, 0]),
    ],
    ids=["preserve", "gene", "custom"],
)
def test_processing_applies_ordering_strategy_and_preserves_values(
    unordered_vcf: Path,
    unordered_catalog: Path,
    ordering: FrequencyOrdering,
    expected_rows: list[int],
) -> None:
    result = calculate_cohort_frequencies(
        unordered_vcf,
        Cohort("test", ("C1-ExC1-xgenV1",)),
        CsvAnnotationSource(unordered_catalog),
        ordering=ordering,
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
    pd.testing.assert_frame_equal(result, expected.iloc[expected_rows])
