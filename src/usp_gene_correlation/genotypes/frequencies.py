"""Validated genotype frequencies per variant."""

from math import isclose

import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from .constants import GENOTYPE_ORDER


class GenotypeFrequencies:
    """Validated frequencies per variant, including the missing-call category."""

    def __init__(self, data: pd.DataFrame) -> None:
        _validate_frequency_columns(data)
        _validate_frequency_values(data)
        self._data = (
            data.loc[:, list(GENOTYPE_ORDER)]
            .copy()
            .rename_axis(index="variant", columns="genotype")
        )

    @property
    def variant_ids(self) -> tuple[str, ...]:
        return tuple(str(variant) for variant in self._data.index)

    def to_dataframe(self) -> pd.DataFrame:
        """Return a copy so callers cannot modify the validated frequencies."""
        return self._data.copy()


def _validate_frequency_columns(data: pd.DataFrame) -> None:
    if data.empty:
        raise ValueError("At least one variant is required for genotype frequencies.")
    if _columns_match_genotype_order(data.columns):
        raise ValueError(f"Frequency columns must be exactly {GENOTYPE_ORDER}.")


def _columns_match_genotype_order(columns: list[str]) -> bool:
    return len(columns) != len(GENOTYPE_ORDER) or set(columns) != set(GENOTYPE_ORDER)


def _validate_frequency_values(data: pd.DataFrame) -> None:
    if not _are_frequencies_numerical(data):
        raise ValueError("Genotype frequencies must be numeric.")
    if _are_frequencies_between_0_and_1(data):
        raise ValueError("Genotype frequencies must be finite values between 0 and 1.")
    if not _frequencies_sum_to_1(data):
        raise ValueError("Genotype frequencies must sum to 1 for each variant.")


def _are_frequencies_numerical(data: pd.DataFrame) -> bool:
    return all(
        is_numeric_dtype(dtype) and not is_bool_dtype(dtype) for dtype in data.dtypes
    )


def _are_frequencies_between_0_and_1(data: pd.DataFrame) -> bool:
    return data.isna().any().any() or not ((data >= 0) & (data <= 1)).all().all()


def _frequencies_sum_to_1(data: pd.DataFrame) -> bool:
    return all(
        isclose(total, 1.0, rel_tol=0, abs_tol=1e-9) for total in data.sum(axis=1)
    )
