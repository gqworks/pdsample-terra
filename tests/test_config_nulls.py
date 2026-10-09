"""Shipped thresholds stay null until a citation is added."""

from pdvar.config import load_config


def test_shipped_config_has_null_penetrance():
    cfg = load_config()
    liability = cfg["liability"]
    assert liability["penetrance"] in ({}, None)
    assert liability["prs_beta"] is None
    assert liability["prs_beta_citation"] is None
    assert cfg["expression"]["min_brain_pext"] is None
    assert cfg["expression"]["unexpressed_weight"] is None
    for locus in cfg["str_loci"].values():
        assert locus["pathogenic_min_repeats"] is None
        assert locus["citation"] is None
    assert "cohort_size" not in cfg
    assert "lookup_max_gnomad_af" in cfg["splicing"]
