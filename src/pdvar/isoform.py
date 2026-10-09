"""Weight LoF and splice variants by brain exon expression.

When the expression table is missing, or a cutoff is unset, the weight is 1
(neutral). A weight below 1 is applied only if the user sets both a cutoff and
``unexpressed_weight``.
"""

from __future__ import annotations

import pandas as pd

from .scoring import NULL_CONSEQUENCES
from .vcf_stream import chrom_bare

LOF_TERMS = set(NULL_CONSEQUENCES) | {"splice_region_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant"}


def _is_weighted(row: pd.Series) -> bool:
    consequence = str(row.get("consequence") or "")
    terms = [part.strip() for part in consequence.replace(",", "&").split("&")]
    if any(term in LOF_TERMS or "splice" in term for term in terms):
        return True
    if str(row.get("dosage_mechanism") or "") == "lof":
        return True
    return bool(row.get("splice_candidate") is True or str(row.get("splice_candidate")) == "True")


def _target(row: pd.Series) -> tuple[str, int, int] | None:
    chrom = chrom_bare(str(row.get("chrom") or ""))
    if not chrom or chrom == "NAN":
        return None
    for start_col, end_col in (
        ("splicevault_exon_start", "splicevault_exon_end"),
        ("pos", "end"),
    ):
        try:
            start = int(row[start_col])
        except (KeyError, TypeError, ValueError):
            continue
        try:
            end = int(row[end_col]) if end_col in row and pd.notna(row[end_col]) else start
        except (TypeError, ValueError):
            end = start
        return chrom, min(start, end), max(start, end)
    return None


def _lookup(expression: pd.DataFrame, chrom: str, start: int, end: int, gene: str) -> tuple[float | None, float | None]:
    gene = gene.upper()
    exp = expression
    if "gene" in exp.columns:
        same = exp["gene"].astype(str).str.upper() == gene
        exp = exp[same] if same.any() else exp.iloc[0:0]
    if exp.empty:
        return None, None
    if {"chrom", "start", "end"}.issubset(exp.columns):
        located = exp.dropna(subset=["chrom", "start", "end"])
        if not located.empty:
            chroms = located["chrom"].map(chrom_bare)
            overlap = located[(chroms == chrom) & (located["start"].astype(int) <= end) & (located["end"].astype(int) >= start)]
            if not overlap.empty:
                exp = overlap
            elif exp["chrom"].notna().any():
                return None, None
    pext = pd.to_numeric(exp["brain_pext"], errors="coerce") if "brain_pext" in exp.columns else pd.Series(dtype=float)
    tpm = (
        pd.to_numeric(exp["gtex_brain_median_tpm"], errors="coerce")
        if "gtex_brain_median_tpm" in exp.columns
        else pd.Series(dtype=float)
    )
    pext_max = None if pext.dropna().empty else float(pext.max())
    tpm_max = None if tpm.dropna().empty else float(tpm.max())
    return pext_max, tpm_max


def _expressed(pext: float | None, tpm: float | None, min_pext, min_tpm) -> bool | None:
    flags = []
    if min_pext is not None and pext is not None:
        flags.append(pext >= float(min_pext))
    if min_tpm is not None and tpm is not None:
        flags.append(tpm >= float(min_tpm))
    if not flags:
        return None
    return any(flags)


def weight_lof(df: pd.DataFrame, expression: pd.DataFrame | None, cfg: dict | None = None) -> pd.DataFrame:
    """Annotate brain expression and set ``lof_weight`` (1 when the resource or cutoff is absent)."""
    cfg = cfg or {}
    df = df.copy()
    df["lof_weight"] = 1.0
    df["lof_weight_source"] = "neutral_missing_resource"
    df["brain_pext"] = pd.NA
    df["gtex_brain_median_tpm"] = pd.NA
    df["brain_expressed"] = pd.NA
    df["weighted_evidence_score"] = df["evidence_score"] if "evidence_score" in df.columns else pd.NA
    if expression is None or expression.empty:
        return df

    min_pext = cfg.get("min_brain_pext")
    min_tpm = cfg.get("min_gtex_tpm")
    downweight = cfg.get("unexpressed_weight")
    for idx, row in df.iterrows():
        target = _target(row)
        if target is None:
            df.at[idx, "lof_weight_source"] = "neutral_no_coordinate"
            continue
        pext, tpm = _lookup(expression, *target, str(row.get("gene") or ""))
        df.at[idx, "brain_pext"] = pd.NA if pext is None else pext
        df.at[idx, "gtex_brain_median_tpm"] = pd.NA if tpm is None else tpm
        expressed = _expressed(pext, tpm, min_pext, min_tpm)
        df.at[idx, "brain_expressed"] = pd.NA if expressed is None else expressed
        if not _is_weighted(row):
            df.at[idx, "lof_weight_source"] = "neutral_not_lof"
            continue
        if expressed is None:
            df.at[idx, "lof_weight_source"] = "annotated_neutral_cutoff_unset"
            continue
        if expressed:
            df.at[idx, "lof_weight_source"] = "brain_expressed"
            continue
        df.at[idx, "lof_weight_source"] = "brain_unexpressed"
        if downweight is not None:
            df.at[idx, "lof_weight"] = float(downweight)
    if "evidence_score" in df.columns:
        df["weighted_evidence_score"] = df["evidence_score"] * df["lof_weight"]
    return df
