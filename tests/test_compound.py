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
