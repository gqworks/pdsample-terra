"""Phenotype-genotype concordance scoring."""

from __future__ import annotations

import pandas as pd


def onset_category(age: pd.Series, pcfg: dict) -> pd.Series:
    cat = pd.Series("unknown", index=age.index, dtype="object")
    cat[age.notna()] = "late"
    cat[age < pcfg["early_onset_age"]] = "early"
    cat[age < pcfg["juvenile_onset_age"]] = "juvenile"
    return cat


def _onset_match(age: float, profile: dict | None) -> object:
    if profile is None or pd.isna(age) or "onset_range" not in profile:
        return pd.NA
    lo, hi = profile["onset_range"]
    return lo <= age <= hi


def _features(text: object, keywords: list[str]) -> list[str]:
    if pd.isna(text) or not keywords:
        return []
    t = str(text).lower()
    return [k for k in keywords if k.lower() in t]


def _ancestry_match(ethnicity: object, profile: dict | None) -> bool:
    if profile is None or pd.isna(ethnicity):
        return False
    e = str(ethnicity).lower()
    return any(a.lower() in e for a in profile.get("enriched_ancestry", []))


def phenotype_concordance(df: pd.DataFrame, pheno: pd.DataFrame, profiles: dict, pcfg: dict) -> pd.DataFrame:
    """Merge patient phenotype and score how well it fits each candidate gene.

    ``phenotype_score`` = +2 onset in expected range, -2 outside, +1 per matched key
    feature, +1 per atypical feature, +1 enriched ancestry.
    """
    df = df.copy()
    df["sample_id"] = df["sample_id"].astype(str)
    df = df.merge(pheno, on="sample_id", how="left", suffixes=("", "_pheno"))
    if "age_onset" in df.columns:
        df["onset_category"] = onset_category(df["age_onset"], pcfg)

    genes = df["gene"].astype(str).str.upper()
    prof = genes.map(lambda g: profiles.get(g))
    age = df["age_onset"] if "age_onset" in df.columns else pd.Series(pd.NA, index=df.index)
    symptoms = df["symptoms"] if "symptoms" in df.columns else pd.Series(pd.NA, index=df.index)
    eth = df["ethnicity"] if "ethnicity" in df.columns else pd.Series(pd.NA, index=df.index)

    df["onset_match"] = pd.array([_onset_match(a, p) for a, p in zip(age, prof)], dtype="boolean")
    df["matched_features"] = [
        ";".join(_features(s, (p or {}).get("key_features", []))) for s, p in zip(symptoms, prof)
    ]
    df["matched_atypical"] = [
        ";".join(_features(s, (p or {}).get("atypical_features", []))) for s, p in zip(symptoms, prof)
    ]
    df["ancestry_enriched"] = [_ancestry_match(e, p) for e, p in zip(eth, prof)]
    df["has_gene_profile"] = prof.notna()

    onset_pts = df["onset_match"].astype(object).map({True: 2, False: -2}).fillna(0).astype(int)
    n_feat = df["matched_features"].map(lambda s: len(s.split(";")) if s else 0)
    n_atyp = df["matched_atypical"].map(lambda s: len(s.split(";")) if s else 0)
    df["phenotype_score"] = onset_pts + n_feat + n_atyp + df["ancestry_enriched"].astype(int)
    return df
