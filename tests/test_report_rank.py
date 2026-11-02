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
