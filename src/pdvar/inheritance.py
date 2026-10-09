"""Mode-of-inheritance logic, including AR second-hit / compound-het detection."""

from __future__ import annotations

import pandas as pd


PANEL_ANNOTATIONS = {
    "moi": "moi",
    "moi_notes": "moi_notes",
    "confidence": "panel_confidence",
    "risk_allele_exception": "risk_allele_exception",
}


def annotate_moi(df: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    """Add panel MOI, notes, confidence and risk-allele flag; genes not on the panel get ``NOT_IN_PANEL``."""
    df = df.copy()
    genes = df["gene"].astype(str).str.upper()
    info = panel.set_index("gene")
    for src, dest in PANEL_ANNOTATIONS.items():
        if src in info.columns:
            df[dest] = genes.map(info[src])
    df["moi"] = df["moi"].fillna("NOT_IN_PANEL") if "moi" in df.columns else "NOT_IN_PANEL"
    return df


def _qualifies(df: pd.DataFrame, min_tier: int) -> pd.Series:
    return df["tier"] <= min_tier


def evaluate_inheritance(df: pd.DataFrame, icfg: dict) -> pd.DataFrame:
    """Assign ``inheritance_status`` per variant, evaluated within each sample x gene.

    Statuses:
      AD_candidate          qualifying variant in an AD gene
      AR_hom                qualifying homozygous variant in an AR gene
      AR_comphet_candidate  >=2 qualifying variants in an AR gene (phase unknown)
      AR_second_hit_splice  the Introme-supported hit of an AR comp-het (partner is AR_comphet_candidate)
      AR_single_het         only one qualifying het in an AR gene (look for missed 2nd hit / CNV)
      XL_hemi               qualifying hemizygous variant in an X-linked gene
      risk_carrier          qualifying variant in a risk-factor gene (e.g. GBA1)
      risk_biallelic        homozygous or >=2 qualifying variants in a risk gene (e.g. GBA1 -> Gaucher)
      MOI_unknown_candidate qualifying variant in a panel gene with inheritance "unknown"
      not_qualifying        below the tier threshold for its MOI, or not on the panel
    """
    df = df.copy()
    df["inheritance_status"] = "not_qualifying"
    df["ar_hit_count"] = 0

    ad_q = _qualifies(df, icfg["ad_min_tier"])
    ar_first = _qualifies(df, icfg["ar_first_hit_min_tier"])
    if "hit_source" in df.columns:
        ar_first &= df["hit_source"] != "introme_lookup"
    ar_second = _qualifies(df, icfg["ar_second_hit_min_tier"])
    if icfg.get("allow_introme_second_hit", True) and "introme_class" in df.columns:
        ar_second |= df["introme_class"] == "high"

    df.loc[(df["moi"] == "AD") & ad_q, "inheritance_status"] = "AD_candidate"
    df.loc[(df["moi"] == "RISK") & ad_q, "inheritance_status"] = "risk_carrier"
    df.loc[(df["moi"] == "XL") & ad_q & (df["zygosity"] == "hemi"), "inheritance_status"] = "XL_hemi"
    df.loc[(df["moi"] == "XL") & ad_q & (df["zygosity"] != "hemi"), "inheritance_status"] = "AD_candidate"
    unknown_q = _qualifies(df, icfg.get("unknown_moi_min_tier", icfg["ad_min_tier"]))
    df.loc[(df["moi"] == "UNKNOWN") & unknown_q, "inheritance_status"] = "MOI_unknown_candidate"

    for (_, _), grp in df[(df["moi"] == "RISK") & ad_q].groupby(["sample_id", "gene"]):
        if len(grp) >= 2 or (grp["zygosity"] == "hom").any():
            df.loc[grp.index, "inheritance_status"] = "risk_biallelic"

    is_ar = df["moi"] == "AR"
    ar_df = df[is_ar & (ar_first | ar_second)]
    for (_, _), grp in ar_df.groupby(["sample_id", "gene"]):
        idx = grp.index
        first = ar_first.loc[idx].to_numpy()
        second = ar_second.loc[idx].to_numpy()
        hom = first & (grp["zygosity"] == "hom").to_numpy()
        n_hits = int(second.sum())
        df.loc[idx, "ar_hit_count"] = n_hits
        if hom.any():
            df.loc[idx[hom], "inheritance_status"] = "AR_hom"
        elif first.any() and n_hits >= 2:
            splice_hit = second & (
                ~first | grp["evidence_reasons"].fillna("").str.contains("introme_high").to_numpy()
            )
            df.loc[idx[second], "inheritance_status"] = "AR_comphet_candidate"
            if splice_hit.any() and (second & ~splice_hit).any():
                df.loc[idx[splice_hit], "inheritance_status"] = "AR_second_hit_splice"
        elif first.any():
            df.loc[idx[first], "inheritance_status"] = "AR_single_het"
    return df
