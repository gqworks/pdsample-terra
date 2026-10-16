"""Phenotype-genotype concordance scoring."""

from __future__ import annotations

import pandas as pd

from pdvar.phenotype_report import looks_like_phenotype_report, parse_phenotype_report

_NUMERIC_FIELDS = {"age_onset", "age_sampling", "age_levodopa_start", "age_diagnosis"}
_FRAME_SKIP = {"comorbidities", "phenotype_source"}


def prepare_phenotype(payload: dict | None) -> dict:
    """Merge a phenotype report with explicit columns. Explicit columns win."""
    payload = dict(payload or {})
    report_text = payload.pop("report_text", None)
    symptoms = payload.get("symptoms")
    if not _filled(report_text) and isinstance(symptoms, str) and looks_like_phenotype_report(symptoms):
        report_text = symptoms
        payload.pop("symptoms", None)
    parsed = parse_phenotype_report(report_text) if _filled(report_text) else {}
    for key, value in payload.items():
        if not _filled(value):
            continue
        parsed[key] = value
    return parsed


def phenotype_for_frame(payload: dict | None) -> dict:
    """Scalar fields that can sit on a candidate row."""
    cleaned = clean_phenotype(prepare_phenotype(payload))
    return {key: value for key, value in cleaned.items() if key not in _FRAME_SKIP}


def _filled(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, dict)):
        return bool(value)
    return True


def clean_phenotype(payload: dict | None) -> dict:
    """Keep non-blank patient fields. Age fields are numeric when they parse.

    Extra keys (biomarkers, family history, a parsed report) are kept so a
    local sample JSON can still feed them. Blank strings are dropped.
    """
    if not payload:
        return {}
    cleaned: dict = {}
    for key, value in payload.items():
        if key in {"sample_id", "sex"} or value is None:
            continue
        if isinstance(value, (list, dict)):
            if value:
                cleaned[key] = value
            continue
        text = str(value).strip()
        if not text:
            continue
        if key in _NUMERIC_FIELDS:
            number = pd.to_numeric(text, errors="coerce")
            if pd.isna(number):
                continue
            number = float(number)
            cleaned[key] = int(number) if number.is_integer() else number
        else:
            cleaned[key] = text
    return cleaned


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


def _text_present(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip()
    return ~text.isin(["", "nan", "None", "<NA>"])


def _clinical_match(provided: bool, profile: bool, positive: bool, onset_out: bool) -> str:
    if not provided:
        return "not_assessed"
    if not profile:
        return "no_profile"
    if positive:
        return "matched"
    if onset_out:
        return "unmatched"
    return "uninformative"


def phenotype_concordance(df: pd.DataFrame, pheno: pd.DataFrame, profiles: dict, pcfg: dict) -> pd.DataFrame:
    """Merge patient phenotype and score how well it fits each candidate gene.

    ``phenotype_score`` = +2 onset in expected range, -2 outside, +1 per matched key
    feature, +1 per atypical feature, +1 enriched ancestry.

    ``clinical_match`` is ``matched`` when onset, a symptom keyword, or ancestry
    agrees with the gene profile. Risk genes without a profile stay ``no_profile``.
    """
    df = df.copy()
    df["sample_id"] = df["sample_id"].astype(str)
    df = df.merge(pheno, on="sample_id", how="left", suffixes=("", "_pheno"))
    if "age_onset" in df.columns:
        df["age_onset"] = pd.to_numeric(df["age_onset"], errors="coerce")
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

    provided = (age.notna() | _text_present(symptoms) | _text_present(eth)).fillna(False)
    ancestry_hit = pd.Series(df["ancestry_enriched"], index=df.index).astype(bool)
    onset_yes = df["onset_match"].eq(True).fillna(False).astype(bool)
    onset_out = df["onset_match"].eq(False).fillna(False).astype(bool)
    positive = (
        onset_yes
        | df["matched_features"].map(bool)
        | df["matched_atypical"].map(bool)
        | ancestry_hit
    ).fillna(False)
    profile = df["has_gene_profile"].fillna(False).astype(bool)
    df["clinical_match"] = [
        _clinical_match(bool(has_input), bool(has_profile), bool(hit), bool(outside))
        for has_input, has_profile, hit, outside in zip(provided, profile, positive, onset_out)
    ]
    return df
