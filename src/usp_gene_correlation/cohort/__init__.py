"""Public API for cohort and CSV loading."""

from .domain import Cohort
from .io import read_cohorts

__all__ = ["Cohort", "read_cohorts"]
