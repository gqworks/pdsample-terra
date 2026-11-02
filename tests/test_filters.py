import pandas as pd

from pdvar import filters

PANEL = pd.DataFrame({"gene": ["PRKN", "LRRK2"]})


def _variants(**overrides):
    base = {
        "sample_id": ["S1"],
        "gene": ["PRKN"],
        "chrom": ["chr6"],
        "pos": [162443325],
        "consequence": ["missense_variant"],
        "zygosity": ["het"],
        "gnomad_af": [0.0001],
        "gnomad_popmax_af": [0.0002],
        "cohort_af": [0.01],
        "depth": [30],
        "gq": [99],
        "allele_balance": [0.5],
        "clinvar_class": ["none"],
        "hit_source": ["small_variant"],
        "variant_class": ["SNV"],
    }
    base.update({k: [v] for k, v in overrides.items()})
    return pd.DataFrame(base)


def _lof_del(**overrides):
    base = {
        "sample_id": ["S1"],
        "gene": ["PRKN"],
        "chrom": ["chr6"],
        "pos": [162386551],
        "end": [162515892],
        "consequence": ["transcript_ablation"],
        "zygosity": ["het"],
        "gnomad_af": [0.0004],
        "gnomad_popmax_af": [0.0004],
        "depth": [None],
        "gq": [None],
        "allele_balance": [None],
        "clinvar_class": ["none"],
        "hit_source": ["sv"],
        "variant_class": ["SV"],
        "svtype": ["DEL"],
        "dosage_event": ["deletion"],
        "dosage_mechanism": ["lof"],
        "scored_dosage": [True],
        "exon_count": [3],
    }
    base.update({k: [v] for k, v in overrides.items()})
    return pd.DataFrame(base)


def test_rare_panel_variant_kept(cfg):
    assert len(filters.apply_filters(_variants(), PANEL, cfg["filters"])) == 1


def test_non_panel_gene_removed(cfg):
    assert filters.apply_filters(_variants(gene="BRCA2"), PANEL, cfg["filters"]).empty


def test_common_variant_removed(cfg):
    assert filters.apply_filters(_variants(gnomad_af=0.2), PANEL, cfg["filters"]).empty


def test_sv_cohort_af_not_hard_filtered(cfg):
    """In-callset SV AF must not drop a rare gnomAD PRKN exon DUP (MH03900-style)."""
    row = _lof_del(
        gnomad_af=2.4e-05,
        gnomad_popmax_af=3e-05,
        sv_panel_af=None,
    )
    # cohort_af mirrors joint VCF INFO/AF (~1.5% in a small PD callset)
    row["cohort_af"] = [0.014925]
    out = filters.apply_filters(row, PANEL, cfg["filters"])
    assert len(out) == 1
    assert out.loc[0, "pass_frequency"]


def test_sv_panel_af_still_hard_filtered(cfg):
    row = _lof_del(gnomad_af=2.4e-05, gnomad_popmax_af=3e-05, sv_panel_af=0.05)
    assert filters.apply_filters(row, PANEL, cfg["filters"]).empty


def test_missing_gnomad_treated_as_rare(cfg):
    out = filters.apply_filters(_variants(gnomad_af=None, gnomad_popmax_af=None), PANEL, cfg["filters"])
    assert len(out) == 1


def test_low_quality_removed(cfg):
    assert filters.apply_filters(_variants(depth=3), PANEL, cfg["filters"]).empty


def test_near_miss_high_impact_depth_rescued(cfg):
    # MH06275-style: frameshift at DP=9 under min_depth=10, within 2x near-miss.
    out = filters.apply_filters(
        _variants(consequence="frameshift_variant", depth=9, zygosity="hom", allele_balance=1.0),
        PANEL,
        cfg["filters"],
    )
    assert len(out) == 1
    assert bool(out.loc[0, "depth_rescued"])
    assert out.loc[0, "depth_rescue_reason"] == "near_miss"


def test_near_miss_does_not_rescue_missense(cfg):
    assert filters.apply_filters(_variants(depth=9), PANEL, cfg["filters"]).empty


def test_near_miss_below_floor_still_removed(cfg):
    assert filters.apply_filters(
        _variants(consequence="frameshift_variant", depth=4, zygosity="hom", allele_balance=1.0),
        PANEL,
        cfg["filters"],
    ).empty


def test_lof_del_rescues_overlapping_low_depth_allele(cfg):
    # Missense would fail near-miss, but overlaps a passing LoF DEL at DP>=5.
    table = pd.concat(
        [
            _variants(depth=6, consequence="missense_variant"),
            _lof_del(),
        ],
        ignore_index=True,
    )
    out = filters.apply_filters(table, PANEL, cfg["filters"])
    rescued = out[out["hit_source"] == "small_variant"]
    assert len(rescued) == 1
    assert bool(rescued.iloc[0]["depth_rescued"])
    assert rescued.iloc[0]["depth_rescue_reason"] == "lof_del"
    assert (out["hit_source"] == "sv").any()


def test_lof_del_does_not_rescue_outside_interval(cfg):
    table = pd.concat(
        [
            _variants(pos=161000000, depth=6, consequence="missense_variant"),
            _lof_del(),
        ],
        ignore_index=True,
    )
    out = filters.apply_filters(table, PANEL, cfg["filters"])
    assert out[out["hit_source"] == "small_variant"].empty


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
