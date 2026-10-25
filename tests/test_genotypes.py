import pandas as pd
import pytest

from usp_gene_correlation.genotypes.domain import GenotypeMatrix, genotype_frequencies


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


def test_matrix_normalizes_calls_and_preserves_axis_order() -> None:
    source = pd.DataFrame(
        {"chr2_2": ["1|0", "0/."], "chr1_1": ["0/0", "1/1"]},
        index=["sample2", "sample1"],
    )
    matrix = GenotypeMatrix(source)
    assert matrix.sample_ids == ("sample2", "sample1")
    assert matrix.variant_ids == ("chr2_2", "chr1_1")
    assert matrix.to_dataframe().loc["sample2", "chr2_2"] == "0/1"
    assert matrix.to_dataframe().loc["sample1", "chr2_2"] == "./."
    assert source.loc["sample2", "chr2_2"] == "1|0"
    assert source.index.name is None
    assert source.columns.name is None
    assert matrix.frequencies().loc["chr2_2"].tolist() == [0, .5, 0, .5]


def test_matrix_isolated_from_input_and_export_mutations() -> None:
    source = pd.DataFrame({"chr1_1": ["0/0"]}, index=["sample"])
    matrix = GenotypeMatrix(source)
    source.iloc[0, 0] = "1/1"
    exported = matrix.to_dataframe()
    exported.iloc[0, 0] = "1/2"
    assert matrix.to_dataframe().iloc[0, 0] == "0/0"
    assert matrix.frequencies().iloc[0].tolist() == [1, 0, 0, 0]


@pytest.mark.parametrize(
    "frame", [pd.DataFrame(), pd.DataFrame({"v": ["1/2"]}), pd.DataFrame({"v": ["1"]})]
)
def test_matrix_rejects_empty_or_unsupported_calls(frame: pd.DataFrame) -> None:
    with pytest.raises(ValueError):
        GenotypeMatrix(frame)
