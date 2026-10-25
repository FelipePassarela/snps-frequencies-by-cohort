"""VCF reading and genotype frequencies."""
from collections.abc import Sequence
from pathlib import Path
from typing import cast
import pandas as pd
import pysam

def genotypes_matrix(
    vcf_path: Path, ids_to_keep: Sequence[str] | None = None
) -> pd.DataFrame:
    """Read VCF and return a DataFrame with samples as rows and variants as columns.

    Create a DataFrame containing the genotypes for each sample (rows) and each variant (columns).
    Genotypes are represented as strings (i.e., "0/0", "0/1", "1/1", "./.").

    Args:
        vcf_path (Path): Path to the VCF file.
        ids_to_keep (Sequence[str]): Sequence of sample IDs to include in the DataFrame. If None, include all samples.

    Returns:
        pd.DataFrame: DataFrame with samples as rows and variants as columns.
    """
    vcf = pysam.VariantFile(str(vcf_path))
    try:
        return _read_genotypes(vcf, ids_to_keep)
    finally:
        vcf.close()


def _read_genotypes(
    vcf: pysam.VariantFile, ids_to_keep: Sequence[str] | None
) -> pd.DataFrame:
    samples: list[str] = list(vcf.header.samples)
    if ids_to_keep is not None:
        ids_set = set(ids_to_keep)  # For O(1) lookups
        samples = [s for s in samples if s in ids_set]

    data: dict[str, list[str]] = {s: [] for s in samples}
    variant_ids: list[str] = []

    for rec in vcf.fetch():
        vid = f"{rec.chrom}_{rec.pos}"
        variant_ids.append(vid)
        for s in samples:
            sample = rec.samples[s]
            gt = cast(tuple[int | None, ...] | None, sample.get("GT"))

            if gt is None:
                gstr = "./."
            else:
                sep = "/"
                # None alleles as '.'
                gstr = sep.join("." if a is None else str(a) for a in gt)
            data[s].append(gstr)

    df = pd.DataFrame(data, index=variant_ids).T
    df.index.name = "sample"
    df.columns.name = "variant"
    return df


def genotype_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """Genotype frequencies for each variant.

    Given a DataFrame with samples as rows and variants as columns,
    return a DataFrame with variants as rows and genotype frequencies as columns.

    Args:
        df (pd.DataFrame): DataFrame with samples as rows and variants as columns.

    Returns:
        pd.DataFrame: DataFrame with variants as rows and genotype frequencies as columns.
    """
    if df.empty:
        raise ValueError("No samples or variants available for frequency calculation.")
    df = df.replace({"1/0": "0/1", "0|1": "0/1", "1|0": "0/1", "0|0": "0/0", "1|1": "1/1"})
    df = df.map(lambda gt: "./." if "." in gt else gt)
    unsupported = set(df.to_numpy().ravel()) - {"0/0", "0/1", "1/1", "./."}
    if unsupported:
        raise ValueError(f"Only diploid biallelic genotypes are supported: {unsupported}")
    counts = df.apply(pd.Series.value_counts).fillna(0).astype(int).T
    totals = counts.sum(axis=1)
    freqs = counts.div(totals, axis=0)

    # Sort columns according to typical genotype order: hom ref > het > hom alt > missing
    freqs = freqs.reindex(columns=["0/0", "0/1", "1/1", "./."], fill_value=0.0)
    freqs.columns.name = "genotype"

    return freqs


