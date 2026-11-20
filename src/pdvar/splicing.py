"""Merge Introme splice predictions and fold them into the evidence score."""

from __future__ import annotations

import pandas as pd

from .scoring import assign_tier


def classify_introme(score: pd.Series, spcfg: dict) -> pd.Series:
    cls = pd.Series("none", index=score.index, dtype="object")
    cls[score >= spcfg["introme_moderate"]] = "moderate"
    cls[score >= spcfg["introme_high"]] = "high"
    return cls


def merge_introme(variants: pd.DataFrame, introme: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in introme.columns if c != "gene"]
    return variants.merge(introme[cols], on="variant_id", how="left")


def apply_introme_evidence(df: pd.DataFrame, spcfg: dict, scfg: dict) -> pd.DataFrame:
    """Add ``introme_class`` and add Introme points to the score for relevant consequences.

    Splice evidence comes from Introme only (VEP SpliceAI is not annotated).
    """
    df = df.copy()
    score = df["introme_score"] if "introme_score" in df.columns else pd.Series(float("nan"), index=df.index)
    df["introme_class"] = classify_introme(score, spcfg)

    consider = set(spcfg["consider_consequences"])
    terms = df["consequence"].fillna("").astype(str).str.split(r"[&,]")
    relevant = terms.map(lambda ts: any(t.strip() in consider for t in ts)).astype(bool)
    reasons = df["evidence_reasons"].fillna("").astype(str) if "evidence_reasons" in df.columns else pd.Series("", index=df.index)
    already = reasons.str.contains("introme_high", regex=False)
    add = (df["introme_class"] == "high") & relevant & ~already

    df["splice_candidate"] = df["introme_class"] != "none"
    df.loc[add, "evidence_score"] += scfg["points"]["introme_high"]
    df.loc[add, "evidence_reasons"] = (reasons.loc[add] + ";introme_high").str.strip(";")
    df["tier"] = assign_tier(df["evidence_score"], scfg["tiers"])
    return df


def ensure_hit_source(df: pd.DataFrame, default: str = "small_variant") -> pd.DataFrame:
    """Fill a missing ``hit_source`` without overwriting SV or callset rows."""
    df = df.copy()
    if "hit_source" not in df.columns:
        df["hit_source"] = default
        return df
    blank = df["hit_source"].isna() | df["hit_source"].astype(str).str.strip().isin(["", "nan", "None", "<NA>"])
    df.loc[blank, "hit_source"] = default
    return df


def summarize_splice_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Record Pangolin and AlphaMissense when the columns already exist. No new points."""
    df = df.copy()
    parts = []
    for column in ("introme_pangolin_gain", "introme_pangolin_loss"):
        if column in df.columns:
            parts.append(pd.to_numeric(df[column], errors="coerce"))
    if parts:
        df["pangolin_max"] = pd.concat(parts, axis=1).max(axis=1)
    else:
        df["pangolin_max"] = pd.NA
    if "alphamissense" not in df.columns:
        df["alphamissense"] = pd.NA
    return df
