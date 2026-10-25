"""Calculate annotated genotype frequencies for a cohort."""

from pathlib import Path

import pandas as pd

from .annotations import SnpAnnotationSource, annotate_frequencies
from .cohort import Cohort
from .genotypes import read_genotypes_from_vcf


def calculate_cohort_frequencies(
    vcf_path: Path, cohort: Cohort, source: SnpAnnotationSource
) -> pd.DataFrame:
    matrix = read_genotypes_from_vcf(vcf_path, ids_to_keep=cohort.sample_ids)
    return annotate_frequencies(matrix.frequencies(), source)
