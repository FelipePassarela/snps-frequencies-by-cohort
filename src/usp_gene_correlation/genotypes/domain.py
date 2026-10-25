"""Genotype normalization, validation and frequencies, independent of VCF I/O."""

import pandas as pd

_GENOTYPE_ORDER = ("0/0", "0/1", "1/1", "./.")
_GENOTYPE_ALIASES = {
    "1/0": "0/1",
    "0|1": "0/1",
    "1|0": "0/1",
    "0|0": "0/0",
    "1|1": "1/1",
}


class GenotypeMatrix:
    """Validated genotype calls with samples as rows and variants as columns."""

    def __init__(self, data: pd.DataFrame) -> None:
        if data.empty:
            raise ValueError(
                "No samples or variants available for frequency calculation."
            )
        normalized = data.map(_normalize_genotype)
        _validate_genotypes(normalized)
        self._data = normalized.rename_axis(index="sample", columns="variant")

    @property
    def sample_ids(self) -> tuple[str, ...]:
        return tuple(str(sample) for sample in self._data.index)

    @property
    def variant_ids(self) -> tuple[str, ...]:
        return tuple(str(variant) for variant in self._data.columns)

    def to_dataframe(self) -> pd.DataFrame:
        """Return a copy so callers cannot modify the validated matrix."""
        return self._data.copy()

    def frequencies(self) -> pd.DataFrame:
        """Calculate frequencies, including missing calls in the denominator."""
        return _calculate_frequencies(self._data)


def genotype_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate frequencies from a DataFrame using the genotype domain model."""
    return GenotypeMatrix(df).frequencies()


def _normalize_genotype(genotype: str) -> str:
    if "." in genotype:
        return "./."
    return _GENOTYPE_ALIASES.get(genotype, genotype)


def _validate_genotypes(genotypes: pd.DataFrame) -> None:
    unsupported = set(genotypes.to_numpy().ravel()) - set(_GENOTYPE_ORDER)
    if unsupported:
        raise ValueError(
            f"Only diploid biallelic genotypes are supported: {unsupported}"
        )


def _calculate_frequencies(genotypes: pd.DataFrame) -> pd.DataFrame:
    counts = genotypes.apply(pd.Series.value_counts).fillna(0).astype(int).T
    totals = counts.sum(axis=1)
    frequencies = counts.div(totals, axis=0)
    frequencies = frequencies.reindex(columns=_GENOTYPE_ORDER, fill_value=0.0)
    frequencies.columns.name = "genotype"
    return frequencies
