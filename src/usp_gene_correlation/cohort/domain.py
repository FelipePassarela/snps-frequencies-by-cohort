"""Cohort domain model, independent of file formats and I/O."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Cohort:
    """A named group of sample identifiers."""

    name: str
    sample_ids: tuple[str, ...]
