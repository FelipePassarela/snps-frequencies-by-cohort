"""Read VCF records and adapt sample calls into a genotype matrix."""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import pandas as pd
import pysam

from .matrix import GenotypeMatrix


@dataclass
class _CollectedGenotypes:
    calls_by_sample: dict[str, list[str]]
    variant_ids: list[str]


def read_genotypes_from_vcf(
    vcf_path: Path, ids_to_keep: Sequence[str] | None = None
) -> GenotypeMatrix:
    """Read samples as rows and variants as columns; None selects all samples."""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        return _read_genotypes(vcf, ids_to_keep)


def _read_genotypes(
    vcf: pysam.VariantFile, ids_to_keep: Sequence[str] | None
) -> GenotypeMatrix:
    samples = _select_samples(list(vcf.header.samples), ids_to_keep)
    collected = _collect_genotypes(vcf, samples)
    return _build_genotype_matrix(collected)


def _select_samples(
    samples: Sequence[str], ids_to_keep: Sequence[str] | None
) -> list[str]:
    """Filter samples while preserving their order in the VCF."""
    if ids_to_keep is None:
        return list(samples)
    selected = set(ids_to_keep)
    return [sample for sample in samples if sample in selected]


def _collect_genotypes(
    vcf: pysam.VariantFile, samples: Sequence[str]
) -> _CollectedGenotypes:
    data: dict[str, list[str]] = {sample: [] for sample in samples}
    variant_ids: list[str] = []

    for record in vcf.fetch():
        variant_ids.append(f"{record.chrom}_{record.pos}")
        for sample in samples:
            data[sample].append(_read_sample_genotype(record, sample))

    return _CollectedGenotypes(calls_by_sample=data, variant_ids=variant_ids)


def _read_sample_genotype(record: pysam.VariantRecord, sample: str) -> str:
    alleles = cast(tuple[int | None, ...] | None, record.samples[sample].get("GT"))
    return _format_genotype(alleles)


def _format_genotype(alleles: tuple[int | None, ...] | None) -> str:
    if alleles is None:
        return "./."
    return "/".join("." if allele is None else str(allele) for allele in alleles)


def _build_genotype_matrix(collected: _CollectedGenotypes) -> GenotypeMatrix:
    result = pd.DataFrame(collected.calls_by_sample, index=collected.variant_ids).T
    result.index.name = "sample"
    result.columns.name = "variant"
    return GenotypeMatrix(result)
