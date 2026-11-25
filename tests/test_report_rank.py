import pandas as pd

from pdvar.report import (
    order_candidate_columns,
    patient_candidates,
    rank_candidates,
    render_html,
    render_markdown,
    write_sample_outputs,
)


def test_hgvs_and_svlen_sit_with_genomic_coordinates():
    frame = pd.DataFrame(
        [
            {
                "consequence": "missense_variant",
                "hgvsp": "ENSP00000000000.1:p.Arg1Cys",
                "rank": 1,
                "variant_id": "1-10-A-G",
                "hgvsc": "ENST00000000000.1:c.1A>G",
                "gene": "PRKN",
                "svlen": pd.NA,
                "chrom": "1",
                "pos": 10,
            }
        ]
    )
    ordered = order_candidate_columns(frame)
    assert list(ordered.columns) == [
        "rank",
        "gene",
        "variant_id",
        "chrom",
        "pos",
        "hgvsc",
        "hgvsp",
        "svlen",
        "consequence",
    ]


def test_followup_reason_is_last_column():
    frame = pd.DataFrame(
        [
            {
                "followup_reason": "single heterozygous hit in recessive gene, second allele not found",
                "rank": 1,
                "gene": "PRKN",
                "inheritance_status": "AR_single_het",
                "variant_id": "6-1-A-G",
                "consequence": "missense_variant",
            }
        ]
    )
    ordered = order_candidate_columns(frame)
    assert list(ordered.columns)[-1] == "followup_reason"
    assert list(ordered.columns)[:4] == ["rank", "gene", "inheritance_status", "variant_id"]


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
            }
        )
    top = patient_candidates(pd.DataFrame(rows), max_per_patient=20)
    ids = set(top["variant_id"])
    assert "missense" in ids
    assert "dup" in ids


def test_followup_single_het_force_kept_after_mendelian_solves():
    """AR single-het follow-up alleles stay in candidates.tsv below solved hits."""
    reason = "single heterozygous hit in recessive gene, second allele not found"
    rows = [
        {
            "sample_id": "S1",
            "gene": "PINK1",
            "variant_id": "solved_hom",
            "consequence": "stop_gained",
            "clinvar_class": "P/LP",
            "inheritance_status": "AR_hom",
            "hit_source": "small_variant",
            "dosage_mechanism": "other",
            "scored_dosage": False,
            "followup_reason": "",
            "evidence_score": 20,
            "tier": 1,
        },
        {
            "sample_id": "S1",
            "gene": "PRKN",
            "variant_id": "single_het",
            "consequence": "missense_variant",
            "clinvar_class": "VUS",
            "inheritance_status": "AR_single_het",
            "hit_source": "small_variant",
            "dosage_mechanism": "other",
            "scored_dosage": False,
            "followup_reason": reason,
            "evidence_score": 5,
            "tier": 3,
        },
    ]
    for i in range(25):
        rows.append(
            {
                "sample_id": "S1",
                "gene": f"GENE{i}",
                "variant_id": f"intron{i}",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "inheritance_status": "AR_hom",
                "hit_source": "small_variant",
                "dosage_mechanism": "other",
                "scored_dosage": False,
                "followup_reason": "",
                "evidence_score": 7,
                "tier": 2,
            }
        )
    top = patient_candidates(pd.DataFrame(rows), max_per_patient=20)
    ids = list(top["variant_id"])
    assert "solved_hom" in ids
    assert "single_het" in ids
    assert ids.index("solved_hom") < ids.index("single_het")
    follow = top[top["variant_id"] == "single_het"].iloc[0]
    assert follow["followup_reason"] == reason


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


def test_ad_gof_intronic_ranks_below_lof_ad_intron_and_str():
    """GoF AD intron noise sits below LoF AD intron demotion and STR expansions."""
    rows = [
        {
            "sample_id": "S1",
            "gene": "ATXN8OS",
            "variant_id": "str",
            "consequence": "short_tandem_repeat_expansion",
            "clinvar_class": "none",
            "inheritance_status": "str_expanded",
            "evidence_score": 8,
            "tier": 2,
        },
        {
            "sample_id": "S1",
            "gene": "GRN",
            "variant_id": "grn_intron",
            "consequence": "intron_variant",
            "clinvar_class": "none",
            "inheritance_status": "AD_candidate",
            "evidence_score": 6,
            "tier": 2,
        },
        {
            "sample_id": "S1",
            "gene": "LRRK2",
            "variant_id": "lrrk2_intron",
            "consequence": "intron_variant",
            "clinvar_class": "none",
            "inheritance_status": "AD_gof_intronic",
            "evidence_score": 20,
            "tier": 1,
        },
    ]
    ranked = rank_candidates(pd.DataFrame(rows))
    assert list(ranked["variant_id"]) == ["str", "grn_intron", "lrrk2_intron"]
    assert ranked.loc[ranked["variant_id"] == "lrrk2_intron", "status_priority"].iloc[0] == 6


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
            }
        )
    ranked = rank_candidates(pd.DataFrame(rows))
    assert ranked.iloc[0]["variant_id"] == "R1441C"
    top = patient_candidates(pd.DataFrame(rows), max_per_patient=20)
    assert top.iloc[0]["variant_id"] == "R1441C"
    assert "R1441C" in set(top["variant_id"])


def test_clinvar_rank_order_plp_then_vus_coding_then_vus_intron():
    """Among same MOI status: P/LP coding > VUS coding > VUS intron > none."""
    rows = [
        {
            "sample_id": "S1",
            "gene": "XPR1",
            "variant_id": "none_missense",
            "consequence": "missense_variant",
            "clinvar_class": "none",
            "inheritance_status": "AD_candidate",
            "evidence_score": 6,
            "tier": 2,
        },
        {
            "sample_id": "S1",
            "gene": "XPR1",
            "variant_id": "vus_intron",
            "consequence": "intron_variant",
            "clinvar_class": "VUS",
            "inheritance_status": "AD_candidate",
            "evidence_score": 3,
            "tier": 3,
        },
        {
            "sample_id": "S1",
            "gene": "XPR1",
            "variant_id": "vus_missense",
            "consequence": "missense_variant",
            "clinvar_class": "VUS",
            "evidence_score": 3,
            "tier": 3,
            "inheritance_status": "AD_candidate",
        },
        {
            "sample_id": "S1",
            "gene": "XPR1",
            "variant_id": "plp_missense",
            "consequence": "missense_variant",
            "clinvar_class": "P/LP",
            "inheritance_status": "AD_candidate",
            "evidence_score": 10,
            "tier": 1,
        },
    ]
    ranked = rank_candidates(pd.DataFrame(rows))
    # Intron AD_candidate gets status_priority 1; coding AD stay at 0.
    assert list(ranked["variant_id"]) == ["plp_missense", "vus_missense", "none_missense", "vus_intron"]


def test_report_copies_unique_qc_json_fields():
    markdown = render_markdown(
        "APGS1222",
        pd.DataFrame(),
        pd.DataFrame(),
        {
            "pass": False,
            "failures": ["FREEMIX 0.05 exceeds 0.02"],
            "sex": {"value": "male", "source": "somalier", "somalier": "male", "provided": "unknown"},
            "contamination": {"freemix": 0.05, "max_freemix": 0.02, "pass": False},
            "coverage": {"mean": 18.2, "min_mean_coverage": 20, "pass": False},
            "reference": {"ok": True, "bam_contigs": 195, "reference_contigs": 195},
            "relatedness": "not_run",
            "modules": {"vep": "called"},
            "novel_gate": "keep-me",
        },
    )
    assert "## QC" in markdown
    assert "pass: False." in markdown
    assert "failures: FREEMIX 0.05 exceeds 0.02." in markdown
    assert "sex: male (source=somalier; somalier=male; provided=unknown)." in markdown
    assert "contamination: FREEMIX 0.05 (max 0.02; pass=False)." in markdown
    assert "coverage: mean 18.2 (min 20; pass=False)." in markdown
    assert '"novel_gate": "keep-me"' in markdown
    assert "## Modules" in markdown


def test_report_omits_not_qualifying_kept_in_tsv(tmp_path):
    candidates = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "PINK1",
                "variant_id": "keep_solve",
                "inheritance_status": "AR_hom",
                "followup_reason": "",
                "rank": 1,
            },
            {
                "sample_id": "S1",
                "gene": "PRKN",
                "variant_id": "rescued_null",
                "inheritance_status": "not_qualifying",
                "followup_reason": "early-onset recessive null",
                "rank": 2,
            },
        ]
    )
    markdown = render_markdown("S1", candidates, pd.DataFrame(), {"modules": {}})
    assert "keep_solve" in markdown
    assert "rescued_null" not in markdown
    html = render_html(markdown, candidates)
    assert "keep_solve" in html
    assert "rescued_null" not in html
    written = write_sample_outputs(candidates, pd.DataFrame(), {}, tmp_path, "S1")
    tsv = written["candidates_tsv"].read_text(encoding="utf-8")
    assert "rescued_null" in tsv
    assert "not_qualifying" in tsv
    html_text = written["report_html"].read_text(encoding="utf-8")
    assert "rescued_null" not in html_text
    assert "rescued_null" not in written["report_md"].read_text(encoding="utf-8")
    assert ">CPRA</th>" in html_text
    assert ">Spliceogenicity</th>" in html_text
    assert ">followup_reason</th>" not in html_text
    assert ">Gene name</th>" in html_text
