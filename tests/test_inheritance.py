import pandas as pd

from pdvar import inheritance


def _run(cfg, rows):
    df = pd.DataFrame(rows, columns=["sample_id", "gene", "moi", "zygosity", "tier", "introme_class", "evidence_reasons"])
    return inheritance.evaluate_inheritance(df, cfg["inheritance"])


def test_ad_candidate(cfg):
    out = _run(cfg, [("S1", "LRRK2", "AD", "het", 1, "none", "clinvar_plp")])
    assert out.loc[0, "inheritance_status"] == "AD_candidate"


def test_ar_homozygous(cfg):
    out = _run(cfg, [("S1", "PRKN", "AR", "hom", 1, "none", "null_variant")])
    assert out.loc[0, "inheritance_status"] == "AR_hom"


def test_ar_single_het(cfg):
    out = _run(cfg, [("S1", "PINK1", "AR", "het", 2, "none", "cadd_high")])
    assert out.loc[0, "inheritance_status"] == "AR_single_het"


def test_ar_comphet(cfg):
    out = _run(cfg, [
        ("S1", "PRKN", "AR", "het", 1, "none", "null_variant"),
        ("S1", "PRKN", "AR", "het", 2, "none", "cadd_high;revel_high"),
    ])
    assert set(out["inheritance_status"]) == {"AR_comphet_candidate"}
    assert (out["ar_hit_count"] == 2).all()


def test_ar_second_hit_from_introme(cfg):
    out = _run(cfg, [
        ("S1", "PRKN", "AR", "het", 1, "none", "null_variant"),
        ("S1", "PRKN", "AR", "het", 4, "high", ""),
    ])
    assert out["inheritance_status"].tolist() == ["AR_comphet_candidate", "AR_second_hit_splice"]


def test_unknown_moi_candidate(cfg):
    out = _run(cfg, [("S1", "DNAJC13", "UNKNOWN", "het", 2, "none", "cadd_high")])
    assert out.loc[0, "inheritance_status"] == "MOI_unknown_candidate"


def test_not_in_panel_never_qualifies(cfg):
    out = _run(cfg, [("S1", "BRCA2", "NOT_IN_PANEL", "het", 1, "none", "null_variant")])
    assert out.loc[0, "inheritance_status"] == "not_qualifying"


def test_risk_gene_biallelic(cfg):
    out = _run(cfg, [
        ("S1", "GBA1", "RISK", "het", 1, "none", "clinvar_plp"),
        ("S1", "GBA1", "RISK", "het", 2, "none", "cadd_high;revel_high"),
        ("S2", "GBA1", "RISK", "het", 1, "none", "clinvar_plp"),
    ])
    assert out["inheritance_status"].tolist() == ["risk_biallelic", "risk_biallelic", "risk_carrier"]


def test_annotate_moi_from_panel():
    panel = pd.DataFrame({"gene": ["GBA1"], "moi": ["RISK"], "confidence": ["high"], "risk_allele_exception": [True]})
    out = inheritance.annotate_moi(pd.DataFrame({"gene": ["GBA1", "BRCA2"]}), panel)
    assert out["moi"].tolist() == ["RISK", "NOT_IN_PANEL"]
    assert out.loc[0, "panel_confidence"] == "high"


def test_gene_alias_mapping():
    from pdvar.normalize import apply_gene_aliases
    assert apply_gene_aliases(pd.Series(["gba", "PRKN"]), {"GBA": "GBA1"}).tolist() == ["GBA1", "PRKN"]


def test_hets_in_different_samples_not_comphet(cfg):
    out = _run(cfg, [
        ("S1", "PRKN", "AR", "het", 1, "none", "null_variant"),
        ("S2", "PRKN", "AR", "het", 1, "none", "null_variant"),
    ])
    assert set(out["inheritance_status"]) == {"AR_single_het"}
