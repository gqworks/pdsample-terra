"""Shipped liability thresholds stay null; STR cutoffs are cited when known."""

from pdvar.config import load_config


def test_shipped_config_has_null_penetrance():
    cfg = load_config()
    liability = cfg["liability"]
    assert liability["penetrance"] in ({}, None)
    assert liability["prs_beta"] is None
    assert liability["prs_beta_citation"] is None
    assert "expression" not in cfg
    assert "cohort_size" not in cfg
    assert "lookup_max_gnomad_af" not in cfg["splicing"]
    assert "lookup_min_score" not in cfg["splicing"]


def test_str_loci_pd_panel_has_cited_cutoffs():
    cfg = load_config()
    loci = cfg["str_loci"]
    # Core EH-catalog PD differentials must carry a numeric floor and a citation.
    required = [
        "ATXN1", "ATXN2", "ATXN3", "CACNA1A", "ATXN7", "ATXN8OS", "ATXN10",
        "PPP2R2B", "TBP", "FMR1", "C9ORF72", "RFC1", "HTT", "AR", "ATN1",
        "NOP56", "JPH3", "GIPC1", "CSTB", "NOTCH2NLC", "FGF14", "FXN",
    ]
    for key in required:
        assert key in loci, key
        assert loci[key]["pathogenic_min_repeats"] is not None, key
        assert loci[key]["citation"], key
    # STRchive now cites numeric floors for SCA31 / OPDM1 (still not in EH catalog).
    for key in ("BEAN1", "LRP12"):
        assert key in loci
        assert loci[key]["pathogenic_min_repeats"] is not None, key
        assert loci[key]["citation"]
