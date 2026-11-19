"""Second-hit rescue, empty-search notes, early-onset keep, and module stubs."""

import pandas as pd

from pdvar import followup, report
from pdvar.compound import NO_EXON_SV, NO_SECOND_SMALL, SINGLE_HET_REASON


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
