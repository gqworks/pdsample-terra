"""Command-line interface: ``pdvar assess``."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from .assess import assess_sample, load_sample_json
from .config import load_config

app = typer.Typer(help="Per-sample Parkinson's disease variant assessment", no_args_is_help=True)


@app.callback()
def main() -> None:
    """Per-sample Parkinson's disease variant assessment.

    A callback keeps ``assess`` as a named subcommand. Typer drops the command
    name when an app has only one command, which makes ``pdvar assess sample.json``
    treat ``assess`` as the sample path.
    """


@app.command()
def assess(
    sample_json: Path = typer.Argument(..., help="Internal per-sample JSON written by the WDL Assess task."),
    config: Optional[Path] = typer.Option(None, "--config", "-c", help="config.yaml"),
    columns: Optional[Path] = typer.Option(None, "--columns", help="columns.yaml"),
    phenotypes: Optional[Path] = typer.Option(None, "--phenotypes", help="gene_phenotypes.yaml"),
    panel: Path = typer.Option(..., "--panel", help="Panel inheritance CSV"),
    gene_bed: Optional[Path] = typer.Option(None, "--gene-bed", help="Gene BED for SV overlap"),
    exon_bed: Optional[Path] = typer.Option(None, "--exon-bed", help="Exon BED for SV overlap"),
    outdir: Path = typer.Option(Path("out"), "--outdir", help="Directory for candidates, reports, and QC JSON"),
) -> None:
    """Rank candidate variants for one sample."""
    cfg = load_config(config, columns, phenotypes)
    sample = load_sample_json(sample_json)
    written = assess_sample(
        sample,
        cfg,
        panel_csv=panel,
        outdir=outdir,
        gene_bed=gene_bed,
        exon_bed=exon_bed,
    )
    for name, path in written.items():
        typer.echo(f"{name}: {path}")


if __name__ == "__main__":
    app()
