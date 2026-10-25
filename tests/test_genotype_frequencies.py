import pandas as pd
import pytest

from usp_gene_correlation.genotypes import GenotypeFrequencies, GenotypeMatrix


def test_matrix_returns_frequency_object() -> None:
    result = GenotypeMatrix(pd.DataFrame({"chr1_1": ["0/0", "./."]})).frequencies()
    assert isinstance(result, GenotypeFrequencies)
    assert result.variant_ids == ("chr1_1",)
    assert result.to_dataframe().iloc[0].tolist() == [0.5, 0, 0, 0.5]


def test_frequencies_preserve_variant_order_and_isolate_data() -> None:
    source = pd.DataFrame(
        {"./.": [0.0, 1.0], "1/1": [0.0, 0.0], "0/1": [0.25, 0.0], "0/0": [0.75, 0.0]},
        index=["chr2_2", "chr1_1"],
    )
    frequencies = GenotypeFrequencies(source)
    assert frequencies.variant_ids == ("chr2_2", "chr1_1")
    assert frequencies.to_dataframe().columns.tolist() == ["0/0", "0/1", "1/1", "./."]
    assert source.index.name is None
    source.iloc[0, 0] = 0.5
    exported = frequencies.to_dataframe()
    exported.iloc[0, 0] = 0.0
    assert frequencies.to_dataframe().loc["chr2_2", "0/0"] == 0.75
    assert frequencies.to_dataframe().loc["chr2_2", "./."] == 0.0


@pytest.mark.parametrize(
    "values",
    [
        [-0.1, 0.6, 0.5, 0.0],
        [1.1, 0.0, 0.0, 0.0],
        [0.2, 0.2, 0.2, 0.2],
        [float("nan"), 0.0, 0.0, 0.0],
        [float("inf"), 0.0, 0.0, 0.0],
        ["1", "0", "0", "0"],
        [True, False, False, False],
    ],
)
def test_invalid_frequency_values(values: list[object]) -> None:
    data = pd.DataFrame([values], columns=["0/0", "0/1", "1/1", "./."])
    with pytest.raises(ValueError):
        GenotypeFrequencies(data)


@pytest.mark.parametrize(
    "columns",
    [
        [],
        ["0/0"],
        ["0/0", "0/1", "1/1", "extra"],
        ["0/0", "0/1", "1/1", "./.", "extra"],
        ["0/0", "0/1", "1/1", "1/1"],
    ],
)
def test_invalid_frequency_columns(columns: list[str]) -> None:
    with pytest.raises(ValueError):
        GenotypeFrequencies(pd.DataFrame([[0.0] * len(columns)], columns=columns))


def test_empty_frequency_table() -> None:
    with pytest.raises(ValueError):
        GenotypeFrequencies(pd.DataFrame(columns=["0/0", "0/1", "1/1", "./."]))


def test_rounding_tolerance() -> None:
    data = pd.DataFrame(
        [[0.1, 0.2, 0.3, 0.40000000001]], columns=["0/0", "0/1", "1/1", "./."]
    )
    GenotypeFrequencies(data)
