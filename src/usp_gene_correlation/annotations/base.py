"""Contract implemented by annotation providers."""

from typing import Protocol

from .models import SnpAnnotation


class SnpAnnotationSource(Protocol):
    def lookup(self, chr_pos: str) -> SnpAnnotation | None: ...
