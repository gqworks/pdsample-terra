"""Rule-based evidence scoring and tiering.

Scores are additive points from ``config.scoring.points``; tiers are assigned from
score cut-offs. This is a triage heuristic, not an ACMG classification.
"""

from __future__ import annotations

import pandas as pd

NULL_CONSEQUENCES = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
}


def is_null_variant(consequence: pd.Series) -> pd.Series:
    terms = consequence.fillna("").astype(str).str.split(r"[&,]")
    return terms.map(lambda ts: any(t.strip() in NULL_CONSEQUENCES for t in ts)).astype(bool)


def _col(df: pd.DataFrame, name: str) -> pd.Series:
    return df[name] if name in df.columns else pd.Series(pd.NA, index=df.index, dtype="Float64")


def _text_col(df: pd.DataFrame, name: str) -> pd.Series:
    if name not in df.columns:
        return pd.Series("", index=df.index, dtype="string")
    return df[name].fillna("").astype(str).str.strip().str.lower()


def evidence_flags(df: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """Return a boolean frame with one column per scoring rule."""
    cadd, revel = _col(df, "cadd"), _col(df, "revel")
    avi = _col(df, "avi_phred")
    eigen = _col(df, "eigen")
    clinvar = df.get("clinvar_class", pd.Series("none", index=df.index))
    sift_pred = _text_col(df, "sift_pred")
    polyphen_pred = _text_col(df, "polyphen_pred")
    muttaster_pred = _text_col(df, "muttaster_pred")
    fathmm_pred = _text_col(df, "fathmm_pred")
    flags = pd.DataFrame(index=df.index)
    flags["clinvar_plp"] = clinvar == "P/LP"
    flags["clinvar_conflicting"] = clinvar == "conflicting"
    flags["clinvar_blb"] = clinvar == "B/LB"
    flags["null_variant"] = is_null_variant(df["consequence"])
    flags["cadd_high"] = (cadd >= scfg["cadd_high"]).fillna(False)
    flags["cadd_moderate"] = ((cadd >= scfg["cadd_moderate"]) & (cadd < scfg["cadd_high"])).fillna(False)
    flags["revel_high"] = (revel >= scfg["revel_high"]).fillna(False)
    # Genome-wide AVI Phred bands (top 0.1% / top 1%), not a clinical cutoff.
    flags["avi_high"] = (avi >= scfg["avi_high"]).fillna(False)
    flags["avi_moderate"] = ((avi >= scfg["avi_moderate"]) & (avi < scfg["avi_high"])).fillna(False)
    flags["sift_deleterious"] = sift_pred.str.contains("deleterious", regex=False) | (sift_pred == "d")
    flags["polyphen_damaging"] = polyphen_pred.isin(
        {"probably_damaging", "possibly_damaging", "d", "p"}
    )
    flags["muttaster_disease"] = muttaster_pred.isin({"a", "d"})
    flags["fathmm_damaging"] = fathmm_pred.isin({"d"})
    flags["eigen_high"] = (eigen >= scfg.get("eigen_high", 10)).fillna(False)
    flags["absent_gnomad"] = _col(df, "gnomad_af").fillna(0) == 0
    return flags.astype(bool)


def assign_tier(score: pd.Series, tiers: dict) -> pd.Series:
    tier = pd.Series(4, index=score.index, dtype="int64")
    tier[score >= tiers["tier3_min_score"]] = 3
    tier[score >= tiers["tier2_min_score"]] = 2
    tier[score >= tiers["tier1_min_score"]] = 1
    return tier


def score_variants(df: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """Add ``evidence_score``, ``evidence_reasons`` and ``tier`` columns."""
    df = df.copy()
    flags = evidence_flags(df, scfg)
    points = scfg["points"]
    rules = [r for r in flags.columns if r in points]
    df["evidence_score"] = sum(flags[r].astype(int) * points[r] for r in rules)
    df["evidence_reasons"] = [
        ";".join(r for r, hit in zip(rules, row) if hit) for row in flags[rules].itertuples(index=False)
    ]
    df["is_null_variant"] = flags["null_variant"]
    df["tier"] = assign_tier(df["evidence_score"], scfg["tiers"])
    return df
