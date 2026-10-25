from pathlib import Path

from usp_gene_correlation.genotypes.io import read_genotypes_from_vcf


def test_vcf_filters_samples_and_normalizes_phasing(vcf: Path) -> None:
    matrix = read_genotypes_from_vcf(vcf, ("C1-ExC1-xgenV1",))
    assert list(matrix.sample_ids) == ["C1-ExC1-xgenV1"]
    assert list(matrix.variant_ids) == ["chr1_1"]
    assert matrix.to_dataframe().iloc[0, 0] == "0/1"
