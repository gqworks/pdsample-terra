"""Default dbnsfp_fields must match academic dbNSFP 5.4a column names."""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COLUMNS = ROOT / "resources" / "dbnsfp" / "dbNSFP5.4a_variant.columns.tsv"
WDL = ROOT / "workflow" / "wdl" / "PdSample.wdl"

# Verified against dbNSFP5.4a_variant.columns.tsv (academic branch).
EXPECTED = [
    "SIFT_score",
    "SIFT_pred",
    "Polyphen2_HDIV_score",
    "Polyphen2_HDIV_pred",
    "Eigen-phred_coding",
    "MutationTaster_score",
    "MutationTaster_pred",
    "fathmm-XF_coding_score",
    "fathmm-XF_coding_pred",
    "CADD_phred",
    "REVEL_score",
]


def test_default_dbnsfp_fields_exist_in_54a_columns_tsv():
    names = {row["name"] for row in csv.DictReader(COLUMNS.open(encoding="utf-8"), delimiter="\t")}
    missing = [name for name in EXPECTED if name not in names]
    assert not missing, f"missing from {COLUMNS.name}: {missing}"
    assert "FATHMM_score" not in names
    assert "FATHMM_pred" not in names


def test_pdsample_default_dbnsfp_fields_string():
    text = WDL.read_text(encoding="utf-8")
    match = re.search(r'String dbnsfp_fields = "([^"]+)"', text)
    assert match, "dbnsfp_fields default not found in PdSample.wdl"
    fields = match.group(1).split(",")
    assert fields == EXPECTED
