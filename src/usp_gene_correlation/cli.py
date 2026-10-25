"""Typer command line interface for cohort genotype frequencies."""

from pathlib import Path
from typing import Annotated, cast

import typer

from .annotations import (
    CsvAnnotationSource,
    SnpAnnotationSource,
    annotate_frequencies,
)
from .cohorts import read_cohorts
from .genotypes import read_genotypes_from_vcf

app = typer.Typer(add_completion=False)


@app.command()
def frequencies(
    ctx: typer.Context,
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
    source = cast(SnpAnnotationSource | None, ctx.obj)
    if source is None:
        source = (
            CsvAnnotationSource(annotations)
            if annotations
            else CsvAnnotationSource.from_default_catalog()
        )

    for cohort in read_cohorts(cohorts):
        typer.echo(f"⚙️ Processing cohort: {cohort.name}")

        matrix = read_genotypes_from_vcf(vcf_path, ids_to_keep=cohort.sample_ids)
        df_frequencies = annotate_frequencies(matrix.frequencies(), source)

        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            cohort_output_path = output_dir / f"{cohort.name}.csv"
            df_frequencies.to_csv(cohort_output_path)
            typer.echo(f"✅ Genotype frequencies saved to {cohort_output_path}")
        else:
            typer.echo(df_frequencies)


def main(
    argv: list[str] | None = None,
    *,
    annotation_source: SnpAnnotationSource | None = None,
) -> None:
    """Run the CLI, preserving programmatic arguments and provider injection."""
    app(args=argv, obj=annotation_source, standalone_mode=argv is None)
