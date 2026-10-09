"""Introme v2 TSV reader. Synthetic rows only."""

from pathlib import Path

import pandas as pd

from pdvar.io import read_table
from pdvar.normalize import normalize_introme

FIXTURES = Path(__file__).parent / "fixtures" / "introme"


def _mapping(cfg):
    return cfg["columns"]["introme"]


def test_raw_v2_tsv_has_scores_and_no_carriers(cfg):
    raw = read_table(FIXTURES / "introme_v2_small.tsv")
    assert "gene" not in raw.columns
    assert "samples_with_alt" not in raw.columns
    variants, carriers = normalize_introme(raw, _mapping(cfg), cfg.get("gene_aliases"))
    assert list(variants["variant_id"]) == ["1-20660000-A-G", "6-161500000-TA-T"]
    row = variants.set_index("variant_id").loc["1-20660000-A-G"]
    assert row["introme_spliceai_max"] == 0.4
    assert row["introme_score"] == 0.85
    assert row["introme_genes"] == ""
    assert row["introme_n_carriers"] == 0
    assert row["introme_n_hom"] == 0
    assert row["introme_cohort_ac"] == 0
    assert carriers.empty


def test_enriched_v2_tsv_accepts_plain_and_literal_genes(cfg):
    raw = read_table(FIXTURES / "introme_v2_enriched.tsv")
    variants, carriers = normalize_introme(raw, _mapping(cfg), cfg.get("gene_aliases"))
    by_id = variants.set_index("variant_id")
    assert by_id.loc["1-20660000-A-G", "introme_genes"] == "PRKN"
    assert by_id.loc["1-20660000-A-G", "introme_n_carriers"] == 1
    assert by_id.loc["1-20660000-A-G", "introme_n_hom"] == 0
    assert by_id.loc["6-161500000-TA-T", "introme_genes"] == "PRKN;PINK1"
    assert by_id.loc["6-161500000-TA-T", "introme_n_hom"] == 1
    assert by_id.loc["1-100-C-T", "introme_genes"] == "PRKN;PINK1"
    assert by_id.loc["1-100-C-T", "introme_n_carriers"] == 0
    assert set(carriers["sample_id"]) == {"SAMPLE_A"}
    hom = carriers[carriers["variant_id"] == "6-161500000-TA-T"]
    assert hom["zygosity"].tolist() == ["hom", "hom"]


def test_plain_alt_is_not_wrapped():
    raw = pd.DataFrame({
        "CHROM": ["chr1"],
        "POS": [10],
        "REF": ["A"],
        "ALT": ["T"],
        "introme_score": [0.9],
    })
    variants, carriers = normalize_introme(raw, {
        "chrom": "CHROM",
        "pos": "POS",
        "ref": "REF",
        "alt": "ALT",
        "introme_score": "introme_score",
    })
    assert variants.loc[0, "variant_id"] == "1-10-A-T"
    assert carriers.empty
