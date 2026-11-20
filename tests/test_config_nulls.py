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
        "BEAN1", "LRP12",
    ]
    for key in required:
        assert key in loci, key
        assert loci[key]["pathogenic_min_repeats"] is not None, key
        assert loci[key]["citation"], key
        assert "not in vendored EH catalog" not in str(loci[key]["citation"])


def test_expansionhunter_catalog_has_twenty_four_loci():
    import json
    from pathlib import Path

    catalog = json.loads(
        (Path(__file__).resolve().parents[1] / "workflow" / "assets" / "expansionhunter_pd.hg38.json").read_text(
            encoding="utf-8"
        )
    )
    ids = [entry["LocusId"] for entry in catalog]
    assert len(ids) == 24
    for locus in ("FGF14", "NOTCH2NLC", "BEAN1", "LRP12", "FXN"):
        assert locus in ids
