"""Second-hit rescue, empty-search notes, early-onset keep, and module stubs."""

import pandas as pd

from pdvar import followup, report
from pdvar.compound import NO_EXON_SV, NO_SECOND_SMALL, SINGLE_HET_REASON, UNPHASED_REASON


def test_rescue_low_depth_null_in_followup_gene():
    scored = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "variant_id": "6-100-C-T",
                "variant_class": "SNV",
                "hit_source": "small_variant",
                "inheritance_status": "AR_single_het",
                "followup_reason": SINGLE_HET_REASON,
                "consequence": "missense_variant",
                "gq": 40,
                "depth": 30,
                "evidence_score": 8,
                "tier": 2,
            }
        ]
    )
    unfiltered = pd.concat(
        [
            scored,
            pd.DataFrame(
                [
                    {
                        "sample_id": "S1",
                        "gene": "PRKN",
                        "variant_id": "6-200-A-AT",
                        "variant_class": "INDEL",
                        "hit_source": "small_variant",
                        "consequence": "frameshift_variant",
                        "clinvar_class": "none",
                        "gq": 25,
                        "depth": 3,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    out = followup.rescue_second_hits(scored, unfiltered, cfg={"filters": {"second_hit_min_depth": 1, "second_hit_min_gq": 20}})
    assert "6-200-A-AT" in set(out["variant_id"].astype(str))
    rescued = out[out["variant_id"] == "6-200-A-AT"].iloc[0]
    assert rescued["inheritance_status"] == "not_qualifying"
    assert "rescued second-hit search" in rescued["followup_reason"]


def test_introme_injection_is_gene_scoped():
    scored = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PINK1",
                "variant_id": "1-10-A-G",
                "followup_reason": SINGLE_HET_REASON,
                "inheritance_status": "AR_single_het",
                "evidence_score": 8,
                "tier": 2,
            }
        ]
    )
    introme = pd.DataFrame(
        [
            {
                "variant_id": "1-99-C-T",
                "introme_genes": "PINK1",
                "introme_score": 0.95,
            },
            {
                "variant_id": "2-99-C-T",
                "introme_genes": "SNCA",
                "introme_score": 0.99,
            },
        ]
    )
    out = followup.rescue_second_hits(scored, scored, introme_variants=introme, cfg={"splicing": {"introme_high": 0.9}}, sample_id="S1")
    ids = set(out["variant_id"].astype(str))
    assert "1-99-C-T" in ids
    assert "2-99-C-T" not in ids


def test_introme_second_hit_keeps_0_9_before_relaxed():
    scored = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PINK1",
                "variant_id": "1-10-A-G",
                "followup_reason": SINGLE_HET_REASON,
                "inheritance_status": "AR_single_het",
                "evidence_score": 8,
                "tier": 2,
            }
        ]
    )
    introme = pd.DataFrame(
        [
            {"variant_id": "1-99-C-T", "introme_genes": "PINK1", "introme_score": 0.95},
            {"variant_id": "1-80-C-T", "introme_genes": "PINK1", "introme_score": 0.75},
        ]
    )
    out = followup.rescue_second_hits(
        scored,
        scored,
        introme_variants=introme,
        cfg={"splicing": {"introme_high": 0.9, "introme_second_hit_relaxed": 0.7}},
        sample_id="S1",
    )
    ids = set(out["variant_id"].astype(str))
    assert "1-99-C-T" in ids
    assert "1-80-C-T" not in ids
    high = out[out["variant_id"] == "1-99-C-T"].iloc[0]
    assert high["introme_class"] == "high"
    assert followup.INTROME_RELAXED_REASON not in str(high["followup_reason"])


def test_introme_second_hit_falls_back_to_0_7_per_gene():
    scored = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PINK1",
                "variant_id": "1-10-A-G",
                "followup_reason": SINGLE_HET_REASON,
                "inheritance_status": "AR_single_het",
                "evidence_score": 8,
                "tier": 2,
            },
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "variant_id": "6-10-A-G",
                "followup_reason": SINGLE_HET_REASON,
                "inheritance_status": "AR_single_het",
                "evidence_score": 8,
                "tier": 2,
            },
        ]
    )
    introme = pd.DataFrame(
        [
            {"variant_id": "6-99-C-T", "introme_genes": "PRKN", "introme_score": 0.95},
            {"variant_id": "1-99-C-T", "introme_genes": "PINK1", "introme_score": 0.75},
            {"variant_id": "1-80-C-T", "introme_genes": "PINK1", "introme_score": 0.65},
            {"variant_id": "2-99-C-T", "introme_genes": "SNCA", "introme_score": 0.75},
        ]
    )
    out = followup.rescue_second_hits(
        scored,
        scored,
        introme_variants=introme,
        cfg={"splicing": {"introme_high": 0.9, "introme_second_hit_relaxed": 0.7}},
        sample_id="S1",
    )
    ids = set(out["variant_id"].astype(str))
    assert "6-99-C-T" in ids
    assert "1-99-C-T" in ids
    assert "1-80-C-T" not in ids
    assert "2-99-C-T" not in ids
    relaxed = out[out["variant_id"] == "1-99-C-T"].iloc[0]
    assert relaxed["introme_class"] == "moderate"
    assert followup.INTROME_RELAXED_REASON in relaxed["followup_reason"]
    assert relaxed["inheritance_status"] == "not_qualifying"


def test_global_introme_high_cutoff_stays_0_9():
    from pdvar.config import load_config

    splicing = load_config()["splicing"]
    assert splicing["introme_high"] == 0.9
    assert splicing["introme_second_hit_relaxed"] == 0.7
    assert splicing["introme_second_hit_relaxed"] < splicing["introme_high"]


def test_empty_search_phrases_and_early_onset_keep():
    df = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "variant_id": "6-100-C-T",
                "variant_class": "SNV",
                "hit_source": "small_variant",
                "inheritance_status": "AR_single_het",
                "followup_reason": SINGLE_HET_REASON,
                "consequence": "missense_variant",
                "clinvar_class": "none",
            }
        ]
    )
    out = followup.annotate_empty_searches(df, gauchian_catalog_negative=False)
    reason = out.iloc[0]["followup_reason"]
    assert NO_SECOND_SMALL in reason
    assert NO_EXON_SV in reason

    nulls = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PINK1",
                "variant_id": "1-50-G-A",
                "inheritance_status": "not_qualifying",
                "consequence": "stop_gained",
                "variant_class": "SNV",
                "hit_source": "small_variant",
                "evidence_score": 8,
                "tier": 4,
            }
        ]
    )
    kept = followup.early_onset_keep(
        nulls,
        age_onset=32,
        classification="unexplained",
        cfg={"phenotype": {"early_onset_age": 50}, "biomarkers": {"gene_classes": {"recessive_eopd": ["PRKN", "PINK1", "PARK7"]}}},
    )
    assert followup.EARLY_ONSET_REASON in kept.iloc[0]["followup_reason"]
    top = report.patient_candidates(kept, max_per_patient=20)
    assert "1-50-G-A" in set(top["variant_id"])

    with_intron_ad = pd.concat(
        [
            nulls,
            pd.DataFrame(
                [
                    {
                        "sample_id": "S1",
                        "gene": "GRN",
                        "variant_id": "17-1-A-G",
                        "inheritance_status": "AD_candidate",
                        "consequence": "intron_variant",
                        "clinvar_class": "none",
                        "mechanism": "lof",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    kept_partial = followup.early_onset_keep(
        with_intron_ad,
        age_onset=32,
        classification="partially_explained",
        cfg={"phenotype": {"early_onset_age": 50}, "biomarkers": {"gene_classes": {"recessive_eopd": ["PRKN", "PINK1", "PARK7"]}}},
    )
    pink1 = kept_partial[kept_partial["variant_id"] == "1-50-G-A"].iloc[0]
    assert followup.EARLY_ONSET_REASON in pink1["followup_reason"]


def test_empty_search_skips_unphased_solved_pairs():
    df = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "variant_id": "6-100-C-T",
                "variant_class": "SNV",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "followup_reason": UNPHASED_REASON,
                "consequence": "missense_variant",
                "clinvar_class": "none",
            },
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "variant_id": "6-200-G-A",
                "variant_class": "SNV",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "followup_reason": UNPHASED_REASON,
                "consequence": "missense_variant",
                "clinvar_class": "none",
            },
        ]
    )
    out = followup.annotate_empty_searches(df)
    assert (out["followup_reason"] == UNPHASED_REASON).all()
    assert NO_EXON_SV not in " ".join(out["followup_reason"].astype(str))


def test_module_stubs_when_required_module_missing():
    stubs = followup.module_stub_rows(
        "S1",
        {"gauchian": "not_assessed", "expansionhunter": "called", "mtdna": "called", "introme": "called", "gatk_sv": "prior"},
        sv_status="prior",
        avi_requested=True,
        pdvar_git_commit="abc123",
    )
    reasons = " ".join(stubs["followup_reason"].astype(str))
    assert "gauchian not assessed" in reasons
    assert "sv_status prior from an older callset" in reasons
    assert "abc123" in reasons
    assert "avi not assessed" in reasons
    assert set(stubs["inheritance_status"]) == {followup.MODULE_STUB_STATUS}
