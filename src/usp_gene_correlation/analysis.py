"""Calculate annotated genotype frequencies for a cohort."""

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .annotations import SnpAnnotationSource, annotate_frequencies
from .cohort import Cohort
from .genotypes import read_genotypes_from_vcf
from .ordering import order_frequencies


def calculate_cohort_frequencies(
    vcf_path: Path,
    cohort: Cohort,
    source: SnpAnnotationSource,
    *,
    gene_order: Sequence[str] | None = None,
) -> pd.DataFrame:
    matrix = read_genotypes_from_vcf(vcf_path, ids_to_keep=cohort.sample_ids)
    result = annotate_frequencies(matrix.frequencies(), source)
    result = order_frequencies(result, gene_order) if gene_order else result
    return result
