from pathlib import Path

import pandas as pd
import pytest

from usp_gene_correlation.ordering import GeneOrder, PreserveOrder, read_gene_order


def test_ordering_preserves_rows_values_and_input(
    annotated_frequencies: pd.DataFrame,
) -> None:
    original = annotated_frequencies.copy(deep=True)

    actual = GeneOrder(("LERFS", "HK1", "ABSENT")).apply(annotated_frequencies)

    pd.testing.assert_frame_equal(actual, original.iloc[[2, 1, 3, 5, 0, 4]])
    pd.testing.assert_frame_equal(annotated_frequencies, original)


def test_ordering_without_matching_genes_preserves_input_order(
    annotated_frequencies: pd.DataFrame,
) -> None:
    pd.testing.assert_frame_equal(
        GeneOrder(("ABSENT",)).apply(annotated_frequencies), annotated_frequencies
    )


def test_preserve_order_returns_original_frequencies(
    annotated_frequencies: pd.DataFrame,
) -> None:
    assert PreserveOrder().apply(annotated_frequencies) is annotated_frequencies


@pytest.mark.parametrize(
    ("contents", "expected_rows"),
    [
        ("name\nLERFS\nHK1\n", [2, 1, 3, 5, 0, 4]),
        ("name\nHK1\nLERFS\n", [1, 3, 5, 2, 0, 4]),
    ],
)
def test_gene_order_csv_defines_result_order(
    tmp_path: Path,
    annotated_frequencies: pd.DataFrame,
    contents: str,
    expected_rows: list[int],
) -> None:
    path = tmp_path / "custom_order.csv"
    path.write_text(contents, encoding="utf-8")

    actual = GeneOrder(read_gene_order(path)).apply(annotated_frequencies)

    pd.testing.assert_frame_equal(actual, annotated_frequencies.iloc[expected_rows])
