import pandas as pd

from pdvar import io, splicing
from pdvar.normalize import normalize_introme


def _parsed(cfg):
    raw = io.load_raw(cfg, "introme")
    return normalize_introme(raw, cfg["columns"]["introme"], cfg.get("gene_aliases"))


def test_parse_introme_variants(cfg):
    variants, _ = _parsed(cfg)
    assert variants["variant_id"].is_unique
    row = variants.set_index("variant_id").loc["6-161500000-TA-T"]
    assert row["introme_genes"] == "ENSG00000286805;PRKN"
    assert row["introme_n_carriers"] == 2  # duplicate rows merged
    common = variants.set_index("variant_id").loc["1-20670000-C-T"]
    assert common["introme_cohort_ac"] == 23
    assert common["introme_n_hom"] == 11
    assert variants.set_index("variant_id").loc["1-20660000-A-G", "introme_spliceai_max"] == 0.4


def test_parse_introme_carriers(cfg):
    _, carriers = _parsed(cfg)
    multi = carriers[carriers["variant_id"] == "6-161500000-TA-T"]
    assert set(multi["gene"]) == {"ENSG00000286805", "PRKN"}
    assert set(multi["sample_id"]) == {"P005", "X1"}
    hom = carriers[(carriers["sample_id"] == "S10")]
    assert hom["zygosity"].tolist() == ["hom"]


def test_find_second_hits_only_for_first_hit_pairs(cfg):
    variants, carriers = _parsed(cfg)
    panel = pd.DataFrame({"gene": ["PINK1", "PRKN"], "moi": ["AR", "AR"]})
    scored = pd.DataFrame({
        "sample_id": ["P003", "P005"],
        "gene": ["PINK1", "PRKN"],
        "variant_id": ["1-20650000-G-A", "6-161360000-C-T"],
        "tier": [2, 4],  # P005 hit does not qualify -> no lookup for P005
    })
    variants = variants.copy()
    variants.loc[variants["variant_id"] == "1-20670000-C-T", "gnomad_af"] = 0.2
    hits = splicing.find_second_hits(scored, carriers, variants, panel, cfg)
    assert hits["variant_id"].tolist() == ["1-20660000-A-G"]
    assert hits["hit_source"].eq("introme_lookup").all()
    assert hits["consequence"].tolist() == ["donor_region"]
    assert hits["introme_class"].tolist() == ["high"]


def test_find_second_hits_af_threshold(cfg):
    variants, carriers = _parsed(cfg)
    panel = pd.DataFrame({"gene": ["PINK1"], "moi": ["AR"]})
    scored = pd.DataFrame({"sample_id": ["P003"], "gene": ["PINK1"], "variant_id": ["x"], "tier": [1]})
    variants = variants.copy()
    variants["gnomad_af"] = 0.2
    hits = splicing.find_second_hits(scored, carriers, variants, panel, cfg)
    assert hits.empty
