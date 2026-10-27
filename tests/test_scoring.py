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


def test_in_silico_predictions_add_points(cfg):
    bare = _score(cfg, consequence="missense_variant", clinvar_class="none", gnomad_af=0.01)
    hit = _score(
        cfg,
        consequence="missense_variant",
        clinvar_class="none",
        gnomad_af=0.01,
        sift_pred="deleterious",
        polyphen_pred="probably_damaging",
        muttaster_pred="D",
        fathmm_pred="D",
        eigen=12.0,
    )
    expected = (
        cfg["scoring"]["points"]["sift_deleterious"]
        + cfg["scoring"]["points"]["polyphen_damaging"]
        + cfg["scoring"]["points"]["muttaster_disease"]
        + cfg["scoring"]["points"]["fathmm_damaging"]
        + cfg["scoring"]["points"]["eigen_high"]
    )
    assert hit["evidence_score"] == bare["evidence_score"] + expected
    for reason in (
        "sift_deleterious",
        "polyphen_damaging",
        "muttaster_disease",
        "fathmm_damaging",
        "eigen_high",
    ):
        assert reason in hit["evidence_reasons"]


def test_avi_phred_bands_add_points(cfg):
    bare = _score(cfg, consequence="missense_variant", clinvar_class="none", gnomad_af=0.0)
    high = _score(cfg, consequence="missense_variant", clinvar_class="none", gnomad_af=0.0, avi_phred=31)
    mid = _score(cfg, consequence="missense_variant", clinvar_class="none", gnomad_af=0.0, avi_phred=21)
    low = _score(cfg, consequence="missense_variant", clinvar_class="none", gnomad_af=0.0, avi_phred=10)
    assert high["evidence_score"] == bare["evidence_score"] + cfg["scoring"]["points"]["avi_high"]
    assert "avi_high" in high["evidence_reasons"]
    assert "avi_moderate" not in high["evidence_reasons"]
    assert mid["evidence_score"] == bare["evidence_score"] + cfg["scoring"]["points"]["avi_moderate"]
    assert "avi_moderate" in mid["evidence_reasons"]
    assert "avi_high" not in mid["evidence_reasons"]
    assert low["evidence_score"] == bare["evidence_score"]
    assert "avi_high" not in low["evidence_reasons"]
    assert "avi_moderate" not in low["evidence_reasons"]


def test_higher_avi_breaks_evidence_ties():
    from pdvar.report import rank_candidates

    frame = pd.DataFrame(
        {
            "sample_id": ["S", "S"],
            "inheritance_status": ["AD_candidate", "AD_candidate"],
            "tier": [2, 2],
            "evidence_score": [6, 6],
            "avi_phred": [12.0, 35.0],
            "variant_id": ["low", "high"],
        }
    )
    ranked = rank_candidates(frame)
    assert ranked.iloc[0]["variant_id"] == "high"


def test_normalize_clinvar():
    s = pd.Series(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_benign", "Conflicting_classifications",
                   "Uncertain_significance", "Benign", None])
    assert normalize_clinvar(s).tolist() == ["P/LP", "P/LP", "other", "conflicting", "VUS", "B/LB", "none"]
