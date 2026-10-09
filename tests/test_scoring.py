import pandas as pd

from pdvar import scoring
from pdvar.normalize import normalize_clinvar


def _score(cfg, **cols):
    df = pd.DataFrame({k: [v] for k, v in cols.items()})
    return scoring.score_variants(df, cfg["scoring"]).iloc[0]


def test_clinvar_pathogenic_is_tier1(cfg):
    row = _score(cfg, consequence="missense_variant", clinvar_class="P/LP", cadd=30.0, revel=0.8, gnomad_af=0.0003)
    assert row["tier"] == 1
    assert "clinvar_plp" in row["evidence_reasons"]


def test_null_variant_with_high_cadd_is_tier1(cfg):
    row = _score(cfg, consequence="frameshift_variant", clinvar_class="none", cadd=35.0, gnomad_af=0.0001)
    assert row["is_null_variant"]
    assert row["tier"] == 1


def test_benign_is_tier4(cfg):
    row = _score(cfg, consequence="synonymous_variant", clinvar_class="B/LB", cadd=2.0, gnomad_af=0.2)
    assert row["tier"] == 4


def test_normalize_clinvar():
    s = pd.Series(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_benign", "Conflicting_classifications",
                   "Uncertain_significance", "Benign", None])
    assert normalize_clinvar(s).tolist() == ["P/LP", "P/LP", "other", "conflicting", "VUS", "B/LB", "none"]
