from pathlib import Path

import pytest

from usp_gene_correlation.cohorts import Cohort, read_cohorts


def test_different_sizes_preserve_order_and_unique_valid_ids(tmp_path: Path) -> None:
    larger = tmp_path / "larger.csv"
    larger.write_text(
        "label;id\n"
        "first;C2-ExC1-xgenV1\n"
        "header;invalid\n"
        "duplicate;C2-ExC1-xgenV1\n"
        "second;C1-ExC1-xgenV1\n"
        "missing;\n"
    )
    smaller = tmp_path / "smaller.csv"
    smaller.write_text("id,label\nC3-ExC1-xgenV1,only\n")
    assert read_cohorts([smaller, larger]) == [
        Cohort("smaller", ("C3-ExC1-xgenV1",)),
        Cohort("larger", ("C2-ExC1-xgenV1", "C1-ExC1-xgenV1")),
    ]


def test_empty_cohort_list() -> None:
    with pytest.raises(ValueError, match="At least one cohort"):
        read_cohorts([])


def test_duplicate_names_rejected_before_reading(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unique stems"):
        read_cohorts([tmp_path / "a/group.csv", tmp_path / "b/group.csv"])


def test_invalid_cohort(tmp_path: Path) -> None:
    path = tmp_path / "invalid.csv"
    path.write_text("name,value\nfoo,bar\n")
    with pytest.raises(ValueError, match="No matching"):
        read_cohorts([path])


def test_read_cohort_filters_header_and_duplicates(cohort: Path) -> None:
    assert read_cohorts([cohort]) == [Cohort("cohort", ("C1-ExC1-xgenV1",))]
