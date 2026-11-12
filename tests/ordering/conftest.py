"""Synthetic annotated frequencies for ordering tests."""

import pandas as pd
import pytest


@pytest.fixture
def annotated_frequencies() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "name": ["OTHER", "HK1", "LERFS", "HK1", "", "HK1"],
            "rsid": ["rs9", "rs2", "rs1", "rs2", "", "rs3"],
            "0/0": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            "0/1": [1.0, 0.8, 0.6, 0.4, 0.2, 0.0],
            "1/1": [0.0] * 6,
            "./.": [0.0] * 6,
        },
        index=pd.Index(
            ["chr1_9", "chr1_2", "chr1_1", "chr1_2", "chr1_8", "chr1_3"],
            name="variant",
        ),
    ).rename_axis(columns="genotype")
