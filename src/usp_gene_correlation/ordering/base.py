"""Contract implemented by frequency ordering strategies."""

from typing import Protocol

import pandas as pd


class FrequencyOrdering(Protocol):
    def apply(self, frequencies: pd.DataFrame) -> pd.DataFrame: ...
