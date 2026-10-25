"""Public API for genotype matrices and frequencies."""

from .matrix import GenotypeMatrix, genotype_frequencies
from .frequencies import GenotypeFrequencies
from .io import read_genotypes_from_vcf

__all__ = ["GenotypeFrequencies", "GenotypeMatrix", "genotype_frequencies", "read_genotypes_from_vcf"]
