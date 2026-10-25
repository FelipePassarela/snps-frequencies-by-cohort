"""Add SNP metadata without changing genotype frequencies."""

from collections.abc import Sequence

import pandas as pd
from loguru import logger

from ..genotypes.frequencies import GenotypeFrequencies
from .base import SnpAnnotationSource
from .models import SnpAnnotation


def annotate_frequencies(
    frequencies: GenotypeFrequencies, source: SnpAnnotationSource
) -> pd.DataFrame:
    """Preserve row order and all frequency values; unknown annotations stay blank."""
    positions = frequencies.variant_ids
    annotations = _lookup_annotations(positions, source)
    result = _with_annotation_columns(frequencies.to_dataframe(), annotations)
    _warn_missing_annotations(positions, annotations)
    return result


def _lookup_annotations(
    positions: Sequence[str], source: SnpAnnotationSource
) -> list[SnpAnnotation | None]:
    return [source.lookup(position) for position in positions]


def _with_annotation_columns(
    frequencies: pd.DataFrame, annotations: Sequence[SnpAnnotation | None]
) -> pd.DataFrame:
    """Insert annotation columns into frequencies."""
    result = frequencies.drop(columns=["name", "rsid"], errors="ignore").copy()
    names = [
        annotation.name if annotation is not None else "" for annotation in annotations
    ]
    rsids = [
        annotation.rsid if annotation is not None else "" for annotation in annotations
    ]
    result.insert(0, "name", names)
    result.insert(1, "rsid", rsids)
    return result


def _warn_missing_annotations(
    positions: Sequence[str], annotations: Sequence[SnpAnnotation | None]
) -> None:
    missing = [
        position
        for position, annotation in zip(positions, annotations)
        if annotation is None
    ]
    if missing:
        logger.warning("No SNP annotation for: {}", ", ".join(missing))
