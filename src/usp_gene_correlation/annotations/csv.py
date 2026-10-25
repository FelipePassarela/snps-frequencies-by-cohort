"""CSV-backed annotation provider."""

import csv
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Self

from .base import SnpAnnotationSource
from .models import SnpAnnotation, normalize_position

_REQUIRED_COLUMNS = ("chr_pos", "name", "rsid")


class CsvAnnotationSource(SnpAnnotationSource):
    """Load a local CSV once; positions must use the VCF genome assembly."""

    def __init__(self, path: Path) -> None:
        self._records = _create_annotations_records(_read_annotations(path))

    def lookup(self, chr_pos: str) -> SnpAnnotation | None:
        return self._records.get(normalize_position(chr_pos))

    @classmethod
    def from_default_catalog(cls) -> Self:
        package_dir = Path(__file__).resolve().parent.parent
        catalog = package_dir / "data" / "snps_annotations.csv"
        return cls(catalog)


def _read_annotations(path: Path) -> list[SnpAnnotation]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        _validate_columns(reader.fieldnames)
        return [_parse_annotation(row) for row in reader]


def _validate_columns(columns: Sequence[str] | None) -> None:
    if not set(_REQUIRED_COLUMNS).issubset(columns or ()):
        raise ValueError("Annotation CSV requires chr_pos,name,rsid columns")


def _parse_annotation(
    row: Mapping[str | None, str | list[str] | None],
) -> SnpAnnotation:
    return SnpAnnotation(
        chr_pos=_required_value(row, "chr_pos"),
        name=_required_value(row, "name"),
        rsid=_required_value(row, "rsid"),
    )


def _required_value(
    row: Mapping[str | None, str | list[str] | None], column: str
) -> str:
    value = row.get(column)
    if not isinstance(value, str):
        raise ValueError("Incomplete annotation CSV row")
    return value.strip()


def _create_annotations_records(
    annotations: Iterable[SnpAnnotation],
) -> dict[str, SnpAnnotation]:
    records: dict[str, SnpAnnotation] = {}
    for annotation in annotations:
        key = normalize_position(annotation.chr_pos)
        _check_duplicated_key(key, records)
        _validate_annotation(annotation, key)
        records[key] = annotation
    return records


def _check_duplicated_key(key: str, records: dict[str, SnpAnnotation]) -> None:
    if key in records:
        raise ValueError(f"Duplicate annotation position: {key}")


def _validate_annotation(annotation: SnpAnnotation, key: str) -> None:
    if not annotation.name.strip() or not re.fullmatch(r"rs[0-9]+", annotation.rsid):
        raise ValueError(f"Invalid annotation for {key}")
