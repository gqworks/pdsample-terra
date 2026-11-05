import pandas as pd

from pdvar.report import patient_candidates, rank_candidates


def test_sv_coding_comphet_outranks_intronic_ar_hom_flood():
    rows = [
        {
            "sample_id": "MH03900",
            "gene": "PRKN",
            "variant_id": "missense",
            "consequence": "missense_variant",
            "clinvar_class": "P/LP",
            "inheritance_status": "AR_comphet_candidate",
            "evidence_score": 21,
            "tier": 1,
            "weighted_evidence_score": 21,
        },
        {
            "sample_id": "MH03900",
            "gene": "PRKN",
            "variant_id": "dup",
            "consequence": "transcript_ablation",
            "clinvar_class": "none",
            "inheritance_status": "AR_comphet_candidate",
            "evidence_score": 8,
            "tier": 2,
            "weighted_evidence_score": 8,
        },
    ]
    for i in range(20):
        rows.append(
            {
                "sample_id": "MH03900",
                "gene": f"GENE{i}",
                "variant_id": f"intron{i}",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AR_hom",
                "evidence_score": 7,
                "tier": 2,
                "weighted_evidence_score": 7,
            }
        )
    df = pd.DataFrame(rows)
    ranked = rank_candidates(df)
    top = patient_candidates(df, max_per_patient=20)
    assert ranked.iloc[0]["variant_id"] == "missense"
    assert ranked.iloc[1]["variant_id"] == "dup"
    assert set(top.head(2)["variant_id"]) == {"missense", "dup"}


def test_plp_and_lof_sv_force_kept_when_buried_by_intron_flood():
    """If pairing never ran, P/LP + LoF SV must still survive the top-N cut."""
    rows = [
        {
            "sample_id": "MH03900",
            "gene": "PRKN",
            "variant_id": "missense",
            "consequence": "missense_variant",
            "clinvar_class": "P/LP",
            "inheritance_status": "AR_single_het",
            "hit_source": "small_variant",
            "dosage_mechanism": "other",
            "scored_dosage": False,
            "evidence_score": 13,
            "tier": 1,
            "weighted_evidence_score": 13,
        },
        {
            "sample_id": "MH03900",
            "gene": "PRKN",
            "variant_id": "dup",
            "consequence": "transcript_ablation",
            "clinvar_class": "none",
            "inheritance_status": "AR_single_het",
            "hit_source": "sv",
            "dosage_mechanism": "lof",
            "scored_dosage": True,
            "evidence_score": 16,
            "tier": 1,
            "weighted_evidence_score": 16,
        },
    ]
    for i in range(25):
        rows.append(
            {
                "sample_id": "MH03900",
                "gene": f"GENE{i}",
                "variant_id": f"intron{i}",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AR_hom",
                "hit_source": "small_variant",
                "dosage_mechanism": "other",
                "scored_dosage": False,
                "evidence_score": 7,
                "tier": 2,
                "weighted_evidence_score": 7,
            }
        )
    top = patient_candidates(pd.DataFrame(rows), max_per_patient=20)
    ids = set(top["variant_id"])
    assert "missense" in ids
    assert "dup" in ids


def test_str_expanded_survives_intronic_ar_flood():
    """ATXN8OS-style expansion must stay in the shortlist vs intron AR priority-0 noise."""
    rows = [
        {
            "sample_id": "MH06257",
            "gene": "ATXN8OS",
            "variant_id": "ATXN8OS_CAG",
            "consequence": "short_tandem_repeat_expansion",
            "clinvar_class": "none",
            "inheritance_status": "str_expanded",
            "hit_source": "str",
            "dosage_mechanism": "other",
            "scored_dosage": False,
            "evidence_score": 8,
            "tier": 2,
            "weighted_evidence_score": 8,
        },
    ]
    for i in range(30):
        rows.append(
            {
                "sample_id": "MH06257",
                "gene": f"GENE{i}",
                "variant_id": f"intron{i}",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AR_hom",
                "hit_source": "small_variant",
                "dosage_mechanism": "other",
                "scored_dosage": False,
                "evidence_score": 7,
                "tier": 2,
                "weighted_evidence_score": 7,
            }
        )
    df = pd.DataFrame(rows)
    # AR_hom (priority 0) still outranks str_expanded (1) before demotion; must-keep saves it.
    top = patient_candidates(df, max_per_patient=20)
    assert "ATXN8OS_CAG" in set(top["variant_id"])
    assert len(top) > 20  # force-kept beyond top-N cut


def test_str_expanded_outranks_demoted_intronic_candidate():
    rows = [
        {
            "sample_id": "MH06257",
            "gene": "ATXN8OS",
            "variant_id": "ATXN8OS_CAG",
            "consequence": "short_tandem_repeat_expansion",
            "clinvar_class": "none",
            "inheritance_status": "str_expanded",
            "evidence_score": 8,
            "tier": 2,
            "weighted_evidence_score": 8,
        },
        {
            "sample_id": "MH06257",
            "gene": "PRKRA",
            "variant_id": "lookup",
            "consequence": "donor_canonical",
            "clinvar_class": "none",
            "inheritance_status": "AR_intronic_candidate",
            "evidence_score": 20,
            "tier": 1,
            "weighted_evidence_score": 20,
        },
    ]
    ranked = rank_candidates(pd.DataFrame(rows))
    assert list(ranked["variant_id"]) == ["ATXN8OS_CAG", "lookup"]


def test_gba1_risk_outranks_intronic_ar_flood():
    """Pathogenic GBA1 risk alleles must not sit under Introme AR_hom intron noise."""
    rows = [
        {
            "sample_id": "MH00246",
            "gene": "GBA1",
            "variant_id": "c.681",
            "consequence": "missense_variant",
            "clinvar_class": "P/LP",
            "inheritance_status": "risk_biallelic",
            "hit_source": "small_variant",
            "evidence_score": 14,
            "tier": 1,
            "weighted_evidence_score": 14,
        },
        {
            "sample_id": "MH00246",
            "gene": "GBA1",
            "variant_id": "c.667",
            "consequence": "missense_variant",
            "clinvar_class": "P/LP",
            "inheritance_status": "risk_biallelic",
            "hit_source": "gauchian",
            "evidence_score": 10,
            "tier": 1,
            "weighted_evidence_score": 10,
        },
    ]
    for i in range(25):
        rows.append(
            {
                "sample_id": "MH00246",
                "gene": f"GENE{i}",
                "variant_id": f"intron{i}",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AR_hom",
                "hit_source": "small_variant",
                "evidence_score": 7,
                "tier": 2,
                "weighted_evidence_score": 7,
            }
        )
    ranked = rank_candidates(pd.DataFrame(rows))
    assert list(ranked.head(2)["variant_id"]) == ["c.681", "c.667"]
    top = patient_candidates(pd.DataFrame(rows), max_per_patient=20)
    assert {"c.681", "c.667"} <= set(top["variant_id"])
    assert list(top.head(2)["variant_id"]) == ["c.681", "c.667"]


def test_introme_ad_candidate_demoted_to_priority_1():
    """Introme intron AD_candidate ranks at priority 1 and is not force-kept."""
    rows = [
        {
            "sample_id": "S1",
            "gene": "LRRK2",
            "variant_id": "R1441C",
            "consequence": "missense_variant",
            "clinvar_class": "P/LP",
            "inheritance_status": "AD_candidate",
            "hit_source": "small_variant",
            "evidence_score": 17,
            "tier": 1,
            "weighted_evidence_score": 17,
        },
        {
            "sample_id": "S1",
            "gene": "PRKN",
            "variant_id": "frameshift",
            "consequence": "frameshift_variant",
            "clinvar_class": "none",
            "inheritance_status": "AR_comphet_candidate",
            "hit_source": "small_variant",
            "evidence_score": 10,
            "tier": 1,
            "weighted_evidence_score": 10,
        },
    ]
    for i in range(25):
        rows.append(
            {
                "sample_id": "S1",
                "gene": f"ADGENE{i}",
                "variant_id": f"introme_ad{i}",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AD_candidate",
                "hit_source": "small_variant",
                "evidence_score": 6,
                "tier": 2,
                "weighted_evidence_score": 6,
            }
        )
    df = pd.DataFrame(rows)
    ranked = rank_candidates(df)
    assert ranked.iloc[0]["variant_id"] == "R1441C"
    assert ranked.iloc[1]["variant_id"] == "frameshift"
    assert (ranked.loc[ranked["variant_id"].str.startswith("introme_ad"), "status_priority"] == 1).all()
    top = patient_candidates(df, max_per_patient=20)
    assert "R1441C" in set(top["variant_id"])
    assert "frameshift" in set(top["variant_id"])
    # Intron AD noise must not force-keep past the top-N cap.
    assert len(top) == 20
    assert sum(vid.startswith("introme_ad") for vid in top["variant_id"]) <= 18


def test_lrrk2_ad_plp_outranks_intron_ar_and_lof_sv_noise():
    """LRRK2 R1441C-style AD P/LP must beat intron AR flood and transcript_ablation SVs."""
    rows = [
        {
            "sample_id": "MH06268",
            "gene": "LRRK2",
            "variant_id": "R1441C",
            "consequence": "missense_variant",
            "clinvar_class": "P/LP",
            "inheritance_status": "AD_candidate",
            "hit_source": "small_variant",
            "evidence_score": 17,
            "tier": 1,
            "weighted_evidence_score": 17,
        },
        {
            "sample_id": "MH06268",
            "gene": "FIG4",
            "variant_id": "fig4_del",
            "consequence": "transcript_ablation",
            "clinvar_class": "none",
            "inheritance_status": "AR_hemizygous",
            "hit_source": "sv",
            "dosage_mechanism": "lof",
            "scored_dosage": True,
            "evidence_score": 18,
            "tier": 1,
            "weighted_evidence_score": 18,
        },
    ]
    for i in range(30):
        rows.append(
            {
                "sample_id": "MH06268",
                "gene": f"GENE{i}",
                "variant_id": f"intron{i}",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AR_comphet_candidate",
                "hit_source": "small_variant",
                "evidence_score": 7,
                "tier": 2,
                "weighted_evidence_score": 7,
            }
        )
    ranked = rank_candidates(pd.DataFrame(rows))
    assert ranked.iloc[0]["variant_id"] == "R1441C"
    top = patient_candidates(pd.DataFrame(rows), max_per_patient=20)
    assert top.iloc[0]["variant_id"] == "R1441C"
    assert "R1441C" in set(top["variant_id"])
