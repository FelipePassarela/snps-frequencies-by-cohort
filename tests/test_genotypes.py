from pathlib import Path

import pandas as pd
import pytest

from usp_gene_correlation.genotypes import genotype_frequencies, genotypes_matrix


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


def test_vcf_filters_samples_and_normalizes_phasing(vcf: Path) -> None:
    matrix = genotypes_matrix(vcf, ("C1-ExC1-xgenV1",))
    assert matrix.index.tolist() == ["C1-ExC1-xgenV1"]
    assert matrix.columns.tolist() == ["chr1_1"]
    assert matrix.iloc[0, 0] == "0/1"
