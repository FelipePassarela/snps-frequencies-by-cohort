"""Calculate annotated genotype frequencies for a cohort."""

from pathlib import Path

import pandas as pd

from .annotations import SnpAnnotationSource, annotate_frequencies
from .cohort import Cohort
from .genotypes import read_genotypes_from_vcf
from .ordering import FrequencyOrdering, PreserveOrder


def calculate_cohort_frequencies(
    vcf_path: Path,
    cohort: Cohort,
    source: SnpAnnotationSource,
    *,
    ordering: FrequencyOrdering = PreserveOrder(),
) -> pd.DataFrame:
    matrix = read_genotypes_from_vcf(vcf_path, ids_to_keep=cohort.sample_ids)
    result = annotate_frequencies(matrix.frequencies(), source)
    return ordering.apply(result)
