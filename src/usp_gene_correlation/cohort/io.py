"""Read cohort definitions from CSV files."""

import re
from pathlib import Path

import pandas as pd

from .domain import Cohort

_SAMPLE_ID_PATTERN = re.compile(r"^C\d+-ExC\d+-xgenV\d+$")


def read_cohorts(paths: list[Path]) -> list[Cohort]:
    """Read one cohort per CSV, preserving file and sample order."""
    _validate_cohort_paths(paths)
    return [_read_cohort(path) for path in paths]


def _validate_cohort_paths(paths: list[Path]) -> None:
    if not paths:
        raise ValueError("At least one cohort is required")
    if len({path.stem for path in paths}) != len(paths):
        raise ValueError("Cohort filenames must have unique stems")


def _read_cohort(path: Path) -> Cohort:
    data = pd.read_csv(path, sep=None, engine="python", header=None)
    return Cohort(name=path.stem, sample_ids=_extract_sample_ids(data))


def _extract_sample_ids(cohort: pd.DataFrame) -> tuple[str, ...]:
    """Use the first column containing sample IDs, discarding headers and duplicates."""
    column = _find_sample_id_column(cohort)
    matches = cohort[column].astype(str).str.match(_SAMPLE_ID_PATTERN)
    return tuple(cohort.loc[matches, column].drop_duplicates().tolist())


def _find_sample_id_column(cohort: pd.DataFrame) -> int:
    # CSVs are read without headers, so their column labels are integer positions.
    for column in cohort.columns:
        if _contains_sample_ids(cohort[column]):
            return int(column)

    raise ValueError(
        f"No matching {_SAMPLE_ID_PATTERN.pattern} columns found in cohort."
    )


def _contains_sample_ids(column: pd.Series) -> bool:
    return bool(column.astype(str).str.match(_SAMPLE_ID_PATTERN).any())
