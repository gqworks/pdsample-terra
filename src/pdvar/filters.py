"""Panel, frequency, quality and consequence filters."""

from __future__ import annotations

import pandas as pd

from .incidental import ttn_truncation_mask


def _le_or_missing(series: pd.Series, threshold: float) -> pd.Series:
    """True where value <= threshold or missing (absent from gnomAD counts as rare)."""
    return series.isna() | (series <= threshold)


def _ge_or_missing(series: pd.Series, threshold: float) -> pd.Series:
    return series.isna() | (series >= threshold)


def panel_mask(df: pd.DataFrame, panel: pd.DataFrame) -> pd.Series:
    return df["gene"].astype(str).str.upper().isin(set(panel["gene"]))


def frequency_mask(df: pd.DataFrame, fcfg: dict) -> pd.Series:
    """gnomAD and, when present, the GATK-SV reference-panel allele frequency.

    There is no in-cohort allele frequency. A sample is one Terra row.
    """
    mask = pd.Series(True, index=df.index)
    for col, key in (
        ("gnomad_af", "max_gnomad_af"),
        ("gnomad_popmax_af", "max_gnomad_popmax_af"),
        ("sv_panel_af", "max_gnomad_af"),
    ):
        if col in df.columns and key in fcfg:
            mask &= _le_or_missing(df[col], fcfg[key])
    return mask


def quality_mask(df: pd.DataFrame, fcfg: dict) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for col, key in (("depth", "min_depth"), ("gq", "min_gq")):
        if col in df.columns:
            mask &= _ge_or_missing(df[col], fcfg[key])
    if "allele_balance" in df.columns:
        het = df.get("zygosity", pd.Series("het", index=df.index)) == "het"
        mask &= ~het | _ge_or_missing(df["allele_balance"], fcfg["min_allele_balance"])
    return mask


def consequence_mask(df: pd.DataFrame, keep: list[str]) -> pd.Series:
    """True if any '&'/','-separated consequence term is in ``keep``."""
    keep_set = set(keep)
    terms = df["consequence"].fillna("").astype(str).str.split(r"[&,]")
    return terms.map(lambda ts: any(t.strip() in keep_set for t in ts)).astype(bool)


def apply_filters(df: pd.DataFrame, panel: pd.DataFrame, fcfg: dict) -> pd.DataFrame:
    """Apply all filters and record which ones each retained variant passed.

    ClinVar P/LP variants are rescued from frequency/consequence filters when
    ``always_keep_clinvar_plp`` is set (they still need to be in the panel and pass QC).
    Panel genes flagged ``risk_allele_exception`` use the relaxed ``risk_exception_*`` AF ceilings.
    """
    df = df.copy()
    df["pass_panel"] = panel_mask(df, panel) if fcfg.get("restrict_to_panel", True) else True
    df["pass_frequency"] = frequency_mask(df, fcfg)
    df["pass_quality"] = quality_mask(df, fcfg)
    df["pass_consequence"] = consequence_mask(df, fcfg["keep_consequences"])

    if "risk_allele_exception" in panel.columns and fcfg.get("risk_exception_max_af") is not None:
        exc_genes = set(panel.loc[panel["risk_allele_exception"].astype(bool), "gene"])
        relaxed = {
            **fcfg,
            "max_gnomad_af": fcfg["risk_exception_max_af"],
            "max_gnomad_popmax_af": fcfg["risk_exception_max_af"],
        }
        exc = df["gene"].astype(str).str.upper().isin(exc_genes) & frequency_mask(df, relaxed)
        df["risk_exception_rescued"] = exc & ~df["pass_frequency"]
        df["pass_frequency"] |= exc

    keep = df["pass_panel"] & df["pass_quality"] & df["pass_frequency"] & df["pass_consequence"]
    sources = fcfg.get("always_keep_sources") or []
    if sources and "hit_source" in df.columns:
        source_keep = df["hit_source"].isin(sources) & df["pass_quality"]
        df["source_rescued"] = source_keep & ~keep
        keep |= source_keep
    if fcfg.get("always_keep_clinvar_plp", True) and "clinvar_class" in df.columns:
        rescued = df["pass_panel"] & df["pass_quality"] & (df["clinvar_class"] == "P/LP")
        df["clinvar_rescued"] = rescued & ~keep
        keep |= rescued
    ttn = ttn_truncation_mask(df)
    df["ttn_truncation"] = ttn
    keep |= ttn
    return df.loc[keep].reset_index(drop=True)
