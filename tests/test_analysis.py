from pathlib import Path

from usp_gene_correlation.analysis import calculate_cohort_frequencies
from usp_gene_correlation.annotations import SnpAnnotation, SnpAnnotationSource
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
