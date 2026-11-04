"""mtDNA-specific evidence scoring and status promotion.

Annotations come from MitImpact 4 INFO tags joined onto the MitochondriaPipeline
VCF. Nuclear CADD/REVEL/AVI rules are not applied here.
"""

from __future__ import annotations

import re

import pandas as pd

# MITOMAP disease status tokens treated as pathogenic / likely pathogenic.
_MITOMAP_PATHOGENIC = re.compile(
    r"(?:Cfrm-\[P\]|Cfrm-\[LP\]|Confirmed.?Pathogenic|Confirmed.?Likely.?Pathogenic)",
    re.I,
)
_APOGEE_PATHOGENIC = {"pathogenic", "likely-pathogenic", "likely_pathogenic"}


def mitomap_pathogenic(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str)
    return text.map(lambda value: bool(_MITOMAP_PATHOGENIC.search(value)))


def apogee_pathogenic(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.lower().str.replace("_", "-", regex=False)
    return text.isin(_APOGEE_PATHOGENIC)


def mitotip_likely_pathogenic(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.upper()
    return text.eq("Q1")


def population_af(df: pd.DataFrame) -> pd.Series:
    """gnomAD population AF for rarity: max of AF_hom and AF_het when either is set.

    Hom alone misses common heteroplasmic alleles. Nuclear rows keep ``gnomad_af``.
    """
    hom = (
        pd.to_numeric(df["gnomad_af"], errors="coerce")
        if "gnomad_af" in df.columns
        else pd.Series(pd.NA, index=df.index, dtype="Float64")
    )
    if "hit_source" not in df.columns or "gnomad_af_het" not in df.columns:
        return hom
    mt = df["hit_source"].fillna("").astype(str).eq("mtdna")
    if not mt.any():
        return hom
    het = pd.to_numeric(df["gnomad_af_het"], errors="coerce")
    combined = pd.concat([hom, het], axis=1).max(axis=1, skipna=True)
    out = hom.copy()
    out.loc[mt] = combined.loc[mt]
    return out


def apply_mtdna_evidence(df: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """Add mtDNA-only evidence points for ``hit_source=mtdna`` rows."""
    if df.empty or "hit_source" not in df.columns:
        return df
    mt = df["hit_source"].fillna("").astype(str).eq("mtdna")
    if not mt.any():
        return df
    df = df.copy()
    points = scfg.get("points") or {}
    if "evidence_reasons" in df.columns:
        reason_list = df["evidence_reasons"].fillna("").astype(str).tolist()
    else:
        reason_list = [""] * len(df)
    score = pd.to_numeric(df.get("evidence_score"), errors="coerce").fillna(0).astype(float)

    def _add(mask: pd.Series, key: str) -> None:
        nonlocal score, reason_list
        bump = float(points.get(key, 0) or 0)
        if bump == 0 or not mask.any():
            return
        hit = (mt & mask).tolist()
        score = score + pd.Series(hit, index=df.index).astype(int) * bump
        for pos, ok in enumerate(hit):
            if not ok:
                continue
            prev = reason_list[pos]
            reason_list[pos] = f"{prev};{key}" if prev else key

    mitomap = (
        mitomap_pathogenic(df["mitomap_status"])
        if "mitomap_status" in df.columns
        else pd.Series(False, index=df.index)
    )
    apogee = (
        apogee_pathogenic(df["apogee2_class"])
        if "apogee2_class" in df.columns
        else pd.Series(False, index=df.index)
    )
    tip = (
        mitotip_likely_pathogenic(df["mitotip_quartile"])
        if "mitotip_quartile" in df.columns
        else pd.Series(False, index=df.index)
    )
    _add(mitomap, "mitomap_pathogenic")
    _add(apogee, "apogee_pathogenic")
    _add(tip, "mitotip_likely_pathogenic")

    df["evidence_score"] = score
    df["evidence_reasons"] = reason_list
    # Re-tier with the same cut-offs used for nuclear alleles.
    from .scoring import assign_tier

    df["tier"] = assign_tier(df["evidence_score"], scfg["tiers"])
    return df


def promote_mtdna_status(df: pd.DataFrame, max_gnomad_af: float | None) -> pd.DataFrame:
    """Upgrade high-evidence mtDNA rows to ``mtdna_pathogenic_candidate``."""
    if df.empty or "hit_source" not in df.columns or "inheritance_status" not in df.columns:
        return df
    df = df.copy()
    mt = df["hit_source"].fillna("").astype(str).eq("mtdna")
    status = df["inheritance_status"].astype(str)
    base = mt & status.isin(["mtdna_variant", "not_qualifying", "nan", ""])
    clin = (
        df["clinvar_class"].fillna("").astype(str).eq("P/LP")
        if "clinvar_class" in df.columns
        else pd.Series(False, index=df.index)
    )
    mitomap = (
        mitomap_pathogenic(df["mitomap_status"])
        if "mitomap_status" in df.columns
        else pd.Series(False, index=df.index)
    )
    tier = pd.to_numeric(df["tier"], errors="coerce") if "tier" in df.columns else pd.Series(pd.NA, index=df.index)
    af = population_af(df)
    if max_gnomad_af is None:
        rare = af.isna() | (af == 0)
    else:
        rare = af.isna() | (af <= float(max_gnomad_af))
    high = clin | mitomap | ((tier <= 2) & rare)
    df.loc[base & high, "inheritance_status"] = "mtdna_pathogenic_candidate"
    return df
