from pathlib import Path

import pytest

from usp_gene_correlation.ordering import read_gene_order


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
