"""Typer command line interface for cohort genotype frequencies."""

from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from .annotations import (
    CsvAnnotationSource,
    SnpAnnotationSource,
)
from .cohorts import read_cohorts
from .analysis import calculate_cohort_frequencies

app = typer.Typer(add_completion=False)


@app.command()
def frequencies(
    vcf_path: Annotated[
        Path,
        typer.Argument(
            exists=True, dir_okay=False, readable=True, help="Path to the VCF file."
        ),
    ],
    cohorts: Annotated[
        list[Path],
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Paths to cohort CSV files.",
        ),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output-dir",
            "--output_dir",
            "-o",
            file_okay=False,
            help="Output directory for frequency CSV files.",
        ),
    ] = None,
    annotations: Annotated[
        Path | None,
        typer.Option(
            "--annotations",
            exists=True,
            dir_okay=False,
            readable=True,
            help="CSV with chr_pos,name,rsid (same genome assembly as VCF).",
        ),
    ] = None,
) -> None:
    """Show genotype frequencies from a VCF file, separated by cohorts."""
    source = _resolve_annotation_source(annotations)
    for cohort in read_cohorts(cohorts):
        typer.echo(f"⚙️ Processing cohort: {cohort.name}")
        result = calculate_cohort_frequencies(vcf_path, cohort, source)
        _output_frequencies(result, cohort.name, output_dir)


def _resolve_annotation_source(catalog_path: Path | None) -> SnpAnnotationSource:
    if catalog_path is not None:
        return CsvAnnotationSource(catalog_path)
    return CsvAnnotationSource.from_default_catalog()


def _output_frequencies(
    frequencies: pd.DataFrame, cohort_name: str, output_dir: Path | None
) -> None:
    if output_dir is None:
        typer.echo(frequencies)
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{cohort_name}.csv"
    frequencies.to_csv(output_path)
    typer.echo(f"✅ Genotype frequencies saved to {output_path}")
