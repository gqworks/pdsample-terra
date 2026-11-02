"""PS phase is preferred for recessive compound hets."""

import pandas as pd

from pdvar.compound import SINGLE_HET_REASON, phase_relation, refine


def _pair(**overrides) -> pd.DataFrame:
    rows = [
        {
            "sample_id": "SAMPLE_A",
            "gene": "PRKN",
            "moi": "AR",
            "variant_id": "6-100-C-T",
            "variant_class": "SNV",
            "chrom": "6",
            "pos": 100,
            "end": 100,
            "ref": "C",
            "alt": "T",
            "genotype": "0/1",
            "phase_set": pd.NA,
            "zygosity": "het",
            "inheritance_status": "AR_comphet_candidate",
        },
        {
            "sample_id": "SAMPLE_A",
            "gene": "PRKN",
            "moi": "AR",
            "variant_id": "6-200-G-A",
            "variant_class": "SNV",
            "chrom": "6",
            "pos": 200,
            "end": 200,
            "ref": "G",
            "alt": "A",
            "genotype": "0/1",
            "phase_set": pd.NA,
            "zygosity": "het",
            "inheritance_status": "AR_comphet_candidate",
        },
    ]
    for row in rows:
        row.update(overrides.get(row["variant_id"], {}))
    return pd.DataFrame(rows)


def test_phase_relation_requires_pipe_and_shared_ps():
    assert phase_relation("0|1", "10", "1|0", "10") == "trans"
    assert phase_relation("0|1", "10", "0|1", "10") == "cis"
    assert phase_relation("0/1", "10", "1/0", "10") == "unphased"
    assert phase_relation("0|1", "10", "1|0", "11") == "unphased"


def test_ps_trans_overrides_unphased_comphet():
    df = _pair(
        **{
            "6-100-C-T": {"genotype": "0|1", "phase_set": "4"},
            "6-200-G-A": {"genotype": "1|0", "phase_set": "4"},
        }
    )
    out, follow = refine(df, {})
    assert set(out["inheritance_status"]) == {"AR_comphet_phased"}
    assert set(out["phase_relation"]) == {"trans"}
    assert set(out["phase_source"]) == {"PS"}
    assert follow.empty


def test_ps_cis_is_a_single_het():
    df = _pair(
        **{
            "6-100-C-T": {"genotype": "0|1", "phase_set": "4"},
            "6-200-G-A": {"genotype": "0|1", "phase_set": "4"},
        }
    )
    out, follow = refine(df, {})
    assert set(out["inheritance_status"]) == {"AR_single_het"}
    assert set(out["phase_source"]) == {"PS"}
    assert follow.iloc[0]["reason"] == SINGLE_HET_REASON


def test_unphased_pair_stays_candidate():
    out, follow = refine(_pair(), {})
    assert set(out["inheritance_status"]) == {"AR_comphet_candidate"}
    assert set(out["phase_source"]) == {""}
    assert follow.empty


def test_hemizygous_overlap_and_homozygous_conflict():
    hemi = pd.DataFrame(
        [
            {
                "sample_id": "SAMPLE_A",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "6-1100-C-T",
                "variant_class": "SNV",
                "chrom": "6",
                "pos": 1100,
                "end": 1100,
                "genotype": "1/1",
                "zygosity": "hom",
                "inheritance_status": "AR_hom",
            },
            {
                "sample_id": "SAMPLE_A",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "sv-del",
                "variant_class": "SV",
                "svtype": "DEL",
                "chrom": "6",
                "pos": 1000,
                "end": 1300,
                "genotype": "0/1",
                "zygosity": "het",
                "copy_number": 1,
                "expected_cn": 2,
                "inheritance_status": "AR_single_het",
            },
        ]
    )
    out, _follow = refine(hemi, {})
    assert set(out["inheritance_status"]) == {"AR_hemizygous"}
    assert set(out["phase_relation"]) == {"hemizygous_overlap"}

    conflict = hemi.copy()
    conflict.loc[conflict["variant_class"] == "SV", "zygosity"] = "hom"
    conflict.loc[conflict["variant_class"] == "SV", "copy_number"] = 0
    conflict.loc[conflict["variant_class"] == "SV", "genotype"] = "1/1"
    out, follow = refine(conflict, {})
    assert set(out["inheritance_status"]) == {"review_conflict"}
    assert "homozygous deletion" in follow.iloc[0]["reason"]


def test_single_het_followup_string():
    df = _pair().iloc[:1].copy()
    df["inheritance_status"] = "AR_single_het"
    _out, follow = refine(df, {})
    assert follow.iloc[0]["reason"] == SINGLE_HET_REASON


def test_prkn_dup_pairs_with_pathogenic_missense_over_intronic_introme():
    """MH03900-style: exon DUP + ClinVar P/LP missense become the AR solve."""
    df = pd.DataFrame(
        [
            {
                "sample_id": "MH03900",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "6-161350208-C-T",
                "variant_class": "SNV",
                "chrom": "6",
                "pos": 161350208,
                "end": 161350208,
                "ref": "C",
                "alt": "T",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "missense_variant",
                "clinvar_class": "P/LP",
                "inheritance_status": "AR_single_het",
                "evidence_score": 21,
                "tier": 1,
            },
            {
                "sample_id": "MH03900",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "dup-exons",
                "variant_class": "SV",
                "svtype": "DUP",
                "hit_source": "sv",
                "chrom": "6",
                "pos": 161959762,
                "end": 162315211,
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "transcript_ablation",
                "dosage_mechanism": "lof",
                "dosage_event": "duplication",
                "scored_dosage": True,
                "clinvar_class": "none",
                "inheritance_status": "AR_single_het",
                "evidence_score": 8,
                "tier": 2,
            },
            {
                "sample_id": "MH03900",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "6-162326073-T-A",
                "variant_class": "SNV",
                "chrom": "6",
                "pos": 162326073,
                "end": 162326073,
                "ref": "T",
                "alt": "A",
                "genotype": "1/1",
                "zygosity": "hom",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AR_hom",
                "evidence_score": 7,
                "tier": 2,
                "introme_class": "high",
            },
        ]
    )
    out, follow = refine(df, {})
    solved = out[out["inheritance_status"] == "AR_comphet_candidate"]
    assert set(solved["variant_id"]) == {"6-161350208-C-T", "dup-exons"}
    assert (out.loc[out["variant_id"] == "6-162326073-T-A", "inheritance_status"] == "AR_intronic_candidate").all()
    assert follow.empty
