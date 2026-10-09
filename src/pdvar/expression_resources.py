"""Build the brain-expression table used by isoform weighting.

GTEx v8 gene-median TPM is a plain download. gnomAD pext is a Hail table
(requester-pays for the v2 base-level matrix) and is accepted here only as a
user-exported TSV. Cutoffs stay null until the user sets them.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

GTEX_V8_GENE_MEDIAN_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)
GNOMAD_PEXT_HELP = "https://gnomad.broadinstitute.org/help/pext"
GNOMAD_PEXT_V2_HAIL = "gs://gnomad-public-requester-pays/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.ht"

BRAIN_TISSUES = ("substantia nigra", "putamen", "caudate", "frontal cortex")

OUTPUT_COLUMNS = [
    "chrom",
    "start",
    "end",
    "gene",
    "gene_id",
    "brain_pext",
    "gtex_brain_median_tpm",
    "expression_resolution",
    "source",
]


def print_urls() -> str:
    return "\n".join(
        [
            f"GTEx v8 gene-median TPM: {GTEX_V8_GENE_MEDIAN_URL}",
            f"gnomAD pext help: {GNOMAD_PEXT_HELP}",
            f"gnomAD v2 base-level pext (Hail, GRCh37, requester-pays, not a direct HTTPS file): {GNOMAD_PEXT_V2_HAIL}",
            "gnomAD v4 pext uses GTEx v10 via gnomAD Hail tables. Export a TSV with columns chrom, start, end, gene, brain_pext.",
            "Brain tissues averaged for the gene-level TPM: substantia nigra, putamen, caudate, frontal cortex.",
        ]
    )


def _brain_columns(columns: list[str]) -> list[str]:
    chosen = []
    for column in columns:
        lowered = column.lower()
        if any(tissue in lowered for tissue in BRAIN_TISSUES):
            chosen.append(column)
    return chosen


def gtex_brain_median(gct_path: str | Path) -> pd.DataFrame:
    """Median TPM across the configured brain tissues. Gene-level, not exon-level."""
    table = pd.read_csv(gct_path, sep="\t", skiprows=2)
    brain = _brain_columns([str(column) for column in table.columns])
    if not brain:
        raise ValueError(f"No substantia nigra, putamen, caudate, or frontal cortex columns in {gct_path}")
    values = table[brain].apply(pd.to_numeric, errors="coerce")
    gene_id = table["Name"] if "Name" in table.columns else pd.Series(pd.NA, index=table.index)
    gene = table["Description"] if "Description" in table.columns else gene_id
    return pd.DataFrame(
        {
            "chrom": pd.NA,
            "start": pd.NA,
            "end": pd.NA,
            "gene": gene.astype(str),
            "gene_id": gene_id.astype(str),
            "brain_pext": pd.NA,
            "gtex_brain_median_tpm": values.median(axis=1),
            "expression_resolution": "gene",
            "source": "gtex_v8_gene_median",
        }
    )


def read_pext(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", comment="#")
    required = {"chrom", "start", "end", "gene", "brain_pext"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"pext TSV is missing {sorted(missing)}")
    out = pd.DataFrame(
        {
            "chrom": table["chrom"].astype(str).str.replace(r"^chr", "", regex=True),
            "start": pd.to_numeric(table["start"], errors="coerce"),
            "end": pd.to_numeric(table["end"], errors="coerce"),
            "gene": table["gene"].astype(str),
            "gene_id": table["gene_id"].astype(str) if "gene_id" in table.columns else pd.NA,
            "brain_pext": pd.to_numeric(table["brain_pext"], errors="coerce"),
            "gtex_brain_median_tpm": pd.NA,
            "expression_resolution": "interval",
            "source": "user_pext_tsv",
        }
    )
    return out


def build_expression_table(gtex: pd.DataFrame | None = None, pext: pd.DataFrame | None = None) -> pd.DataFrame:
    frames = [frame for frame in (gtex, pext) if frame is not None and not frame.empty]
    if not frames:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    return pd.concat(frames, ignore_index=True)[OUTPUT_COLUMNS]
