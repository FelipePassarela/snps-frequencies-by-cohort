"""Strategies for preserving VCF order or sorting by gene priority."""

from dataclasses import dataclass

import pandas as pd

from .base import FrequencyOrdering


class PreserveOrder(FrequencyOrdering):
    """Leave the frequency rows in their original order."""

    def apply(self, frequencies: pd.DataFrame) -> pd.DataFrame:
        return frequencies


@dataclass(frozen=True)
class GeneOrder(FrequencyOrdering):
    """Keep unlisted genes last and preserve input order within each gene."""

    genes: tuple[str, ...]

    def apply(self, frequencies: pd.DataFrame) -> pd.DataFrame:
        priorities = {name: rank for rank, name in enumerate(self.genes)}
        return frequencies.sort_values(
            by="name",
            key=lambda names: names.map(priorities).fillna(len(self.genes)),
            kind="stable",
        )
