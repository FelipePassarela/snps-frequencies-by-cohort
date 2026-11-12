"""Read dedicated gene orders and sort annotated frequencies without filtering."""

from collections.abc import Sequence
from pathlib import Path

import pandas as pd


def read_gene_order(path: Path) -> tuple[str, ...]:
    """Read a CSV with a single name column, rejecting empty or repeated genes."""
    genes = _read_gene_names(path)
    _validate_gene_names(genes)
    return tuple(genes)


def _read_gene_names(path: Path) -> pd.Series:
    try:
        data = pd.read_csv(
            path,
            dtype=str,
            keep_default_na=False,
            skip_blank_lines=False,
            encoding="utf-8-sig",
        )
    except pd.errors.EmptyDataError as error:
        raise ValueError("Gene order CSV requires a single name column.") from error
    except pd.errors.ParserError as error:
        raise ValueError(f"Invalid gene order CSV: {error}") from error
    if data.columns.tolist() != ["name"]:
        raise ValueError("Gene order CSV requires a single name column.")
    # Extra fields can be interpreted by pandas as an implicit row index.
    if not isinstance(data.index, pd.RangeIndex):
        raise ValueError("Gene order CSV requires one nonempty name per row.")

    return data["name"].str.strip()


def _validate_gene_names(genes: pd.Series) -> None:
    if genes.empty:
        raise ValueError("Gene order CSV must contain at least one gene.")
    if genes.eq("").any():
        raise ValueError("Gene order CSV requires one nonempty name per row.")
    duplicates = genes[genes.duplicated()]
    if not duplicates.empty:
        raise ValueError(f"Duplicate gene in order CSV: {duplicates.iloc[0]}")


def order_frequencies(
    frequencies: pd.DataFrame, gene_order: Sequence[str]
) -> pd.DataFrame:
    """Keep unlisted genes last and preserve input order within each gene."""
    priorities = {name: rank for rank, name in enumerate(gene_order)}
    return frequencies.sort_values(
        by="name",
        key=lambda names: names.map(priorities).fillna(len(gene_order)),
        kind="stable",
    )
