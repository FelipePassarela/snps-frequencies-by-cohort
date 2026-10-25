"""Public annotation API, independent of the underlying provider."""

from .base import SnpAnnotationSource
from .csv import CsvAnnotationSource
from .enrichment import annotate_frequencies
from .models import SnpAnnotation, normalize_position

__all__ = [
    "SnpAnnotationSource",
    "CsvAnnotationSource",
    "SnpAnnotation",
    "annotate_frequencies",
    "normalize_position",
]
