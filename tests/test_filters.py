import pandas as pd

from pdvar import filters

PANEL = pd.DataFrame({"gene": ["PRKN", "LRRK2"]})


def _variants(**overrides):
    base = {
        "gene": ["PRKN"],
        "consequence": ["missense_variant"],
        "zygosity": ["het"],
        "gnomad_af": [0.0001],
        "gnomad_popmax_af": [0.0002],
        "cohort_af": [0.01],
        "depth": [30],
        "gq": [99],
        "allele_balance": [0.5],
        "clinvar_class": ["none"],
    }
    base.update({k: [v] for k, v in overrides.items()})
    return pd.DataFrame(base)


def test_rare_panel_variant_kept(cfg):
    assert len(filters.apply_filters(_variants(), PANEL, cfg["filters"])) == 1


def test_non_panel_gene_removed(cfg):
    assert filters.apply_filters(_variants(gene="BRCA2"), PANEL, cfg["filters"]).empty


def test_common_variant_removed(cfg):
    assert filters.apply_filters(_variants(gnomad_af=0.2), PANEL, cfg["filters"]).empty


def test_missing_gnomad_treated_as_rare(cfg):
    out = filters.apply_filters(_variants(gnomad_af=None, gnomad_popmax_af=None), PANEL, cfg["filters"])
    assert len(out) == 1


def test_low_quality_removed(cfg):
    assert filters.apply_filters(_variants(depth=3), PANEL, cfg["filters"]).empty


def test_clinvar_plp_rescued_from_frequency(cfg):
    out = filters.apply_filters(_variants(gnomad_af=0.02, clinvar_class="P/LP"), PANEL, cfg["filters"])
    assert len(out) == 1
    assert out.loc[0, "clinvar_rescued"]


def test_risk_allele_exception_relaxes_af(cfg):
    panel = pd.DataFrame({"gene": ["PRKN", "GBA1"], "risk_allele_exception": [False, True]})
    common = {"gnomad_af": 0.02, "gnomad_popmax_af": 0.03}
    assert len(filters.apply_filters(_variants(gene="GBA1", **common), panel, cfg["filters"])) == 1
    assert filters.apply_filters(_variants(gene="PRKN", **common), panel, cfg["filters"]).empty
    assert filters.apply_filters(_variants(gene="GBA1", gnomad_af=0.2), panel, cfg["filters"]).empty


def test_consequence_multi_term():
    df = pd.DataFrame({"consequence": ["splice_region_variant&intron_variant", "downstream_gene_variant"]})
    assert filters.consequence_mask(df, ["splice_region_variant"]).tolist() == [True, False]
