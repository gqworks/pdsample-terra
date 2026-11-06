"""Recessive second hits from overlapping deletions and coding partners."""

import pandas as pd

from pdvar.compound import SINGLE_HET_REASON, refine


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
            "zygosity": "het",
            "inheritance_status": "AR_comphet_candidate",
        },
    ]
    for row in rows:
        row.update(overrides.get(row["variant_id"], {}))
    return pd.DataFrame(rows)


def test_unphased_pair_stays_candidate():
    out, follow = refine(_pair())
    assert set(out["inheritance_status"]) == {"AR_comphet_candidate"}
    assert "phase_source" not in out.columns
    assert "phase_relation" not in out.columns
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
    out, _follow = refine(hemi)
    assert set(out["inheritance_status"]) == {"AR_hemizygous"}

    conflict = hemi.copy()
    conflict.loc[conflict["variant_class"] == "SV", "zygosity"] = "hom"
    conflict.loc[conflict["variant_class"] == "SV", "copy_number"] = 0
    conflict.loc[conflict["variant_class"] == "SV", "genotype"] = "1/1"
    out, follow = refine(conflict)
    assert set(out["inheritance_status"]) == {"review_conflict"}
    assert "homozygous deletion" in follow.iloc[0]["reason"]


def test_single_het_followup_string():
    df = _pair().iloc[:1].copy()
    df["inheritance_status"] = "AR_single_het"
    _out, follow = refine(df)
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
    out, follow = refine(df)
    solved = out[out["inheritance_status"] == "AR_comphet_candidate"]
    assert set(solved["variant_id"]) == {"6-161350208-C-T", "dup-exons"}
    assert (out.loc[out["variant_id"] == "6-162326073-T-A", "inheritance_status"] == "AR_intronic_candidate").all()
    assert follow.empty


def test_intron_only_ar_gene_demoted_to_intronic_candidate():
    """Genes whose only AR 'solves' are intron / splice-region noise drop to support status."""
    df = pd.DataFrame(
        [
            {
                "sample_id": "MH06257",
                "gene": "PRKRA",
                "moi": "AR",
                "variant_id": "2-178450414-A-AT",
                "variant_class": "INDEL",
                "chrom": "2",
                "pos": 178450414,
                "end": 178450414,
                "ref": "A",
                "alt": "AT",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "splice_region_variant&intron_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "evidence_score": 7,
                "tier": 2,
            },
            {
                "sample_id": "MH06257",
                "gene": "PRKRA",
                "moi": "AR",
                "variant_id": "2-178450500-C-T",
                "variant_class": "SNV",
                "chrom": "2",
                "pos": 178450500,
                "end": 178450500,
                "ref": "C",
                "alt": "T",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "evidence_score": 7,
                "tier": 2,
            },
            {
                "sample_id": "MH06257",
                "gene": "NOISE1",
                "moi": "AR",
                "variant_id": "1-100-A-G",
                "variant_class": "SNV",
                "chrom": "1",
                "pos": 100,
                "end": 100,
                "ref": "A",
                "alt": "G",
                "genotype": "1/1",
                "zygosity": "hom",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_hom",
                "evidence_score": 7,
                "tier": 2,
            },
        ]
    )
    out, follow = refine(df)
    assert (out["inheritance_status"] == "AR_intronic_candidate").all()
    assert follow.empty


def test_demoted_intron_hom_does_not_leave_false_comphet():
    """Single coding het must not stay AR_comphet after intron AR_hom demotion."""
    df = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "frameshift",
                "variant_class": "INDEL",
                "chrom": "6",
                "pos": 100,
                "end": 100,
                "ref": "AT",
                "alt": "A",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "frameshift_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "evidence_score": 10,
                "tier": 1,
            },
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "intron_hom",
                "variant_class": "SNV",
                "chrom": "6",
                "pos": 200,
                "end": 200,
                "ref": "A",
                "alt": "G",
                "genotype": "1/1",
                "zygosity": "hom",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_hom",
                "evidence_score": 6,
                "tier": 2,
                "introme_class": "high",
            },
        ]
    )
    out, follow = refine(df)
    by_id = out.set_index("variant_id")["inheritance_status"]
    assert by_id["frameshift"] == "AR_single_het"
    assert by_id["intron_hom"] == "AR_intronic_candidate"
    assert not follow.empty


def test_introme_ar_hom_demoted_when_coding_partners_exist():
    """MONOPDAUS_000050-style: coding comp-het unmasks Introme intron AR_hom."""
    df = pd.DataFrame(
        [
            {
                "sample_id": "MONOPDAUS_000050",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "6-162443378-CCT-C",
                "variant_class": "INDEL",
                "chrom": "6",
                "pos": 162443378,
                "end": 162443378,
                "ref": "CCT",
                "alt": "C",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "frameshift_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "evidence_score": 10,
                "tier": 1,
            },
            {
                "sample_id": "MONOPDAUS_000050",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "6-161569384-G-A",
                "variant_class": "SNV",
                "chrom": "6",
                "pos": 161569384,
                "end": 161569384,
                "ref": "G",
                "alt": "A",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "missense_variant",
                "clinvar_class": "VUS",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "evidence_score": 7,
                "tier": 2,
            },
            {
                "sample_id": "MONOPDAUS_000050",
                "gene": "PRKN",
                "moi": "AR",
                "variant_id": "6-162201629-A-G",
                "variant_class": "SNV",
                "chrom": "6",
                "pos": 162201629,
                "end": 162201629,
                "ref": "A",
                "alt": "G",
                "genotype": "1/1",
                "zygosity": "hom",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_hom",
                "evidence_score": 6,
                "tier": 2,
                "introme_class": "high",
            },
        ]
    )
    out, follow = refine(df)
    by_id = out.set_index("variant_id")["inheritance_status"]
    assert by_id["6-162443378-CCT-C"] == "AR_comphet_candidate"
    assert by_id["6-161569384-G-A"] == "AR_comphet_candidate"
    assert by_id["6-162201629-A-G"] == "AR_intronic_candidate"
    assert follow.empty


def test_coding_plus_intron_keeps_coding_solved():
    """Coding allele keeps solved status; gene-level demotion does not fire when strong exists."""
    df = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PARK7",
                "moi": "AR",
                "variant_id": "missense",
                "variant_class": "SNV",
                "chrom": "1",
                "pos": 10,
                "end": 10,
                "ref": "A",
                "alt": "G",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "missense_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "evidence_score": 10,
                "tier": 1,
            },
            {
                "sample_id": "S1",
                "gene": "PARK7",
                "moi": "AR",
                "variant_id": "intron",
                "variant_class": "SNV",
                "chrom": "1",
                "pos": 20,
                "end": 20,
                "ref": "C",
                "alt": "T",
                "genotype": "0/1",
                "zygosity": "het",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "hit_source": "small_variant",
                "inheritance_status": "AR_comphet_candidate",
                "evidence_score": 7,
                "tier": 2,
            },
        ]
    )
    out, _follow = refine(df)
    by_id = out.set_index("variant_id")["inheritance_status"]
    assert by_id["missense"] == "AR_comphet_candidate"
    # Strong allele present → intron partner is not gene-level demoted.
    assert by_id["intron"] == "AR_comphet_candidate"
