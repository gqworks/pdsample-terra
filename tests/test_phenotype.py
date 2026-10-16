"""Phenotype concordance for risk and Mendelian candidates."""

import pandas as pd

from pdvar.config import load_config
from pdvar.phenotype import clean_phenotype, phenotype_concordance
from pdvar.report import render_markdown


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
    assert out.loc["SMPD1", "clinical_match"] == "no_profile"


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
