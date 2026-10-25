"""SNP annotation value object and position normalization."""

from dataclasses import dataclass


@dataclass(frozen=True)
class SnpAnnotation:
    chr_pos: str
    name: str
    rsid: str


def normalize_position(value: str) -> str:
    """Return a chromosome-position key with a chr prefix and a positive position."""
    chromosome, position = _parse_position(value.strip().removeprefix("chr"))
    return f"chr{chromosome}_{position}"


def _parse_position(value: str) -> tuple[str, int]:
    chromosome, separator, position = value.rpartition("_")
    if not separator or not chromosome or not position.isdigit():
        raise ValueError(f"Invalid chromosome position: {value!r}")

    numeric_position = int(position)
    if numeric_position < 1:
        raise ValueError(f"Invalid chromosome position: {value!r}")

    return chromosome, numeric_position
