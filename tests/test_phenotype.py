"""Phenotype concordance for risk and Mendelian candidates."""

import csv
from pathlib import Path

import pandas as pd

from pdvar.biomarkers import apply_biomarker_priors
from pdvar.config import load_config
from pdvar.phenotype import clean_phenotype, phenotype_concordance, phenotype_for_frame, prepare_phenotype
from pdvar.phenotype_report import parse_phenotype_report
from pdvar.report import render_markdown

GP2_REPORT = """\
GP2ID: APGS_002617 | STUDYCODE_ID: parkinson03684
Sex: Male | Age at sample collection: 58 | PD status: Yes
Levodopa: Yes | Age at levodopa start: 57 | Age at diagnosis: NA | Diagnosis type: Neurologist | Age at onset: 56

Exposures and comorbidities

High clinical relevance:
- High blood pressure (age: 50.0)
- High cholesterol (age: 50.0)

Lower relevance:
- constipation (age: 48)
- cold sores (age: 12)
- back problems (age: 35)
- Hearing loss (age: 58.0)
"""


def _score(rows, pheno_row):
    cfg = load_config()
    return phenotype_concordance(
        pd.DataFrame(rows),
        pd.DataFrame([pheno_row]),
        cfg["gene_phenotypes"],
        cfg["phenotype"],
    )


def test_string_age_matches_gba1_and_flags_genes_without_a_profile():
    out = _score(
        [
            {"sample_id": "S", "gene": "GBA1", "inheritance_status": "risk_carrier"},
            {"sample_id": "S", "gene": "PRKN", "inheritance_status": "AR_single_het"},
            {"sample_id": "S", "gene": "SMPD1", "inheritance_status": "risk_carrier"},
            {"sample_id": "S", "gene": "NOTAPANELGENE", "inheritance_status": "risk_carrier"},
        ],
        {
            "sample_id": "S",
            "age_onset": "55",
            "symptoms": "cognitive impairment and dementia",
            "ethnicity": "Ashkenazi Jewish",
        },
    ).set_index("gene")
    gba = out.loc["GBA1"]
    assert bool(gba["onset_match"]) is True
    assert "dementia" in gba["matched_features"]
    assert "cognitive impairment" in gba["matched_features"]
    assert bool(gba["ancestry_enriched"]) is True
    assert gba["phenotype_score"] == 5
    assert gba["clinical_match"] == "matched"
    assert out.loc["PRKN", "clinical_match"] == "unmatched"
    assert out.loc["PRKN", "phenotype_score"] == -2
    assert bool(out.loc["SMPD1", "ancestry_enriched"]) is True
    assert out.loc["SMPD1", "clinical_match"] == "matched"
    assert out.loc["NOTAPANELGENE", "clinical_match"] == "no_profile"


def test_blank_phenotype_is_not_assessed_and_sparse_text_is_uninformative():
    blank = _score(
        [{"sample_id": "S", "gene": "GBA1"}],
        {"sample_id": "S", "sex": "female", "age_onset": "  ", "symptoms": ""},
    )
    assert blank.iloc[0]["clinical_match"] == "not_assessed"
    sparse = _score(
        [{"sample_id": "S", "gene": "GBA1"}],
        {"sample_id": "S", "symptoms": "fatigue"},
    )
    assert sparse.iloc[0]["clinical_match"] == "uninformative"
    assert sparse.iloc[0]["phenotype_score"] == 0


def test_clean_phenotype_parses_age_and_drops_blanks():
    assert clean_phenotype({"age_onset": "55", "symptoms": " tremor ", "ethnicity": ""}) == {
        "age_onset": 55,
        "symptoms": "tremor",
    }
    assert clean_phenotype(None) == {}
    assert clean_phenotype({"age_onset": "unknown", "saa": "positive"}) == {"saa": "positive"}


def test_report_names_a_missing_phenotype_and_lists_risk_rows():
    missing = render_markdown("S", pd.DataFrame(), pd.DataFrame(), {"modules": {}, "phenotype": {"provided": False}})
    assert "clinical concordance was not assessed" in missing
    present = render_markdown(
        "S",
        pd.DataFrame(
            [
                {
                    "rank": 1,
                    "gene": "GBA1",
                    "moi": "RISK",
                    "inheritance_status": "risk_carrier",
                    "variant_id": "1-10-A-G",
                    "clinical_match": "matched",
                    "onset_match": True,
                    "matched_features": "dementia",
                    "matched_atypical": "",
                    "ancestry_enriched": True,
                    "phenotype_score": 4,
                }
            ]
        ),
        pd.DataFrame(),
        {"modules": {}, "phenotype": {"provided": True, "age_onset": 55, "symptoms": "dementia", "ethnicity": "Ashkenazi Jewish"}},
    )
    section = present.split("## Phenotype concordance", 1)[1].split("## Ranked candidates", 1)[0]
    assert "Age at onset: 55." in section
    assert "GBA1" in section
    assert "matched" in section


def test_gp2_phenotype_report_sets_onset_match_and_keeps_comorbidities():
    parsed = parse_phenotype_report(GP2_REPORT)
    assert parsed["gp2_id"] == "APGS_002617"
    assert parsed["study_id"] == "parkinson03684"
    assert parsed["reported_sex"] == "male"
    assert parsed["pd_status"] == "Yes"
    assert parsed["levodopa"] == "Yes"
    assert parsed["diagnosis_type"] == "Neurologist"
    assert "age_diagnosis" not in parsed
    high = [item for item in parsed["comorbidities"] if item["relevance"] == "high"]
    lower = [item for item in parsed["comorbidities"] if item["relevance"] == "lower"]
    assert [item["condition"] for item in high] == ["High blood pressure", "High cholesterol"]
    assert high[0]["age"] == 50
    assert [item["condition"] for item in lower] == ["constipation", "cold sores", "back problems", "Hearing loss"]
    assert lower[0]["age"] == 48

    pasted = clean_phenotype(prepare_phenotype({"symptoms": GP2_REPORT}))
    assert pasted["age_onset"] == 56
    assert pasted["age_sampling"] == 58
    assert pasted["age_levodopa_start"] == 57
    assert "symptoms" not in pasted

    fields = phenotype_for_frame({"report_text": GP2_REPORT, "age_onset": "40"})
    assert fields["age_onset"] == 40
    assert "comorbidities" not in fields

    fields = phenotype_for_frame({"report_text": GP2_REPORT})
    out = _score(
        [
            {"sample_id": "S", "gene": "GBA1", "inheritance_status": "risk_carrier"},
            {"sample_id": "S", "gene": "PRKN", "inheritance_status": "AR_single_het"},
        ],
        {"sample_id": "S", **fields},
    ).set_index("gene")
    assert out.loc["GBA1", "clinical_match"] == "matched"
    assert bool(out.loc["GBA1", "onset_match"]) is True
    assert out.loc["GBA1", "matched_features"] == ""
    assert out.loc["GBA1", "phenotype_score"] == 2
    assert out.loc["PRKN", "clinical_match"] == "unmatched"
    assert out.loc["PRKN", "matched_features"] == ""

    noted = apply_biomarker_priors(
        pd.DataFrame([{"sample_id": "S", "gene": "PRKN", "levodopa": fields["levodopa"], "age_onset": fields["age_onset"]}]),
        load_config()["biomarkers"],
        early_onset_age=50,
    )
    assert "levodopa_response:flag_only" in noted.iloc[0]["biomarker_notes"]
    assert noted.iloc[0]["biomarker_prior_adjusted"] == False  # noqa: E712

    info = clean_phenotype(prepare_phenotype({"report_text": GP2_REPORT}))
    report = render_markdown("S", pd.DataFrame(), pd.DataFrame(), {"modules": {}, "phenotype": {"provided": True, **info}})
    section = report.split("## Phenotype concordance", 1)[1]
    assert "APGS_002617" in section
    assert "Age at diagnosis: blank." in section
    assert "Comorbidities (high): High blood pressure (age 50); High cholesterol (age 50)." in section
    assert "Comorbidities (lower): constipation (age 48); cold sores (age 12); back problems (age 35); Hearing loss (age 58)." in section


REPO = Path(__file__).resolve().parents[1]
PROFILE_FIELDS = {"moi", "onset_range", "key_features", "atypical_features", "enriched_ancestry"}


def test_every_panel_gene_has_a_phenotype_profile():
    profiles = load_config()["gene_phenotypes"]
    with open(REPO / "resources/panel/monopdaus-gene-inheritance.csv", newline="") as handle:
        symbols = [row["Symbol"] for row in csv.DictReader(handle)]
    expected = {("GBA1" if s == "GBA" else s).upper() for s in symbols}
    assert expected - set(profiles) == set()
    for gene, profile in profiles.items():
        assert gene == gene.upper(), gene
        assert set(profile) <= PROFILE_FIELDS, gene
        assert {"moi", "key_features", "atypical_features", "enriched_ancestry"} <= set(profile), gene
        if "onset_range" in profile:
            lo, hi = profile["onset_range"]
            assert 0 <= lo <= hi <= 100, gene
        for field in ("key_features", "atypical_features", "enriched_ancestry"):
            assert all(isinstance(k, str) and k.strip() for k in profile[field]), (gene, field)


def test_profile_without_onset_range_is_never_unmatched_on_age():
    out = _score(
        [{"sample_id": "S", "gene": "DNAJC13"}],
        {"sample_id": "S", "age_onset": "90", "symptoms": "tremor"},
    )
    assert out.iloc[0]["clinical_match"] == "uninformative"
    assert out.iloc[0]["phenotype_score"] == 0
