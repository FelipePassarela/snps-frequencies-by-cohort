"""Public API for genotype matrices and frequencies."""

from .domain import GenotypeMatrix, genotype_frequencies
from .io import read_genotypes_from_vcf

__all__ = ["GenotypeMatrix", "genotype_frequencies", "read_genotypes_from_vcf"]
