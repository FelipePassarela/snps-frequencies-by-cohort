from pathlib import Path

import pandas as pd
import pytest

from usp_gene_correlation.ordering import GeneOrder, PreserveOrder, read_gene_order


def test_read_gene_order_handles_bom_and_trims_names(tmp_path: Path) -> None:
    path = tmp_path / "order.csv"
    path.write_text("name\n HK1 \nLERFS\n", encoding="utf-8-sig")

    assert read_gene_order(path) == ("HK1", "LERFS")


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("", "single name column"),
        ("gene\nHK1\n", "single name column"),
        ("name,group\nHK1,longa\n", "single name column"),
        ("name\n", "at least one gene"),
        ("name\nHK1\n \n", "one nonempty name"),
        ("name\nHK1,LERFS\n", "one nonempty name"),
        ("name\nHK1\n HK1 \n", "Duplicate gene"),
        ('name\n"HK1', "Invalid gene order CSV"),
    ],
)
def test_read_gene_order_rejects_invalid_files(
    tmp_path: Path, contents: str, message: str
) -> None:
    path = tmp_path / "invalid_order.csv"
    path.write_text(contents)

    with pytest.raises(ValueError, match=message):
        read_gene_order(path)


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
    ("filename", "expected"),
    [
        (
            "longa_gene_order.csv",
            (
                "LERFS",
                "HK1",
                "CPNE9",
                "LIPA",
                "HINT3",
                "IMPA2",
                "CNGB1",
                "SYNE1",
                "ANTXR2",
                "JAKMIP3",
                "PNKD",
                "CMYA5",
                "PIK3R5",
                "NKAIN4",
                "VASH1",
                "TRPV4",
                "CENPJ",
                "LAMB4",
                "GSDMD",
            ),
        ),
        (
            "grave_gene_order.csv",
            (
                "LINC01118",
                "SLC23A2",
                "CCSER1",
                "TRIM6",
                "ARRDC2",
                "SNED1",
                "PODN",
                "SLC4A2",
                "RNF222",
                "ITPR1",
                "TMEM11",
                "COL19A1",
                "C3",
                "GAD2",
                "LASP1",
                "PTPN5",
                "PLA2G4D",
                "RAI14",
                "NPC1",
                "SCN9A",
                "ODAD1",
                "COLEC11",
                "PPARGC1A",
                "NRP2",
                "MUC16",
            ),
        ),
    ],
)
def test_provided_gene_order_files_match_requested_sequences(
    filename: str, expected: tuple[str, ...]
) -> None:
    path = Path(__file__).resolve().parents[1] / "data" / "genotypes" / filename

    assert read_gene_order(path) == expected
