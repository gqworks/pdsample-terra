"""Unit tests for academic dbNSFP 5.4a per-chrom → VEP-ready prep."""

from __future__ import annotations

import gzip
import importlib.util
import shutil
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "prepare_dbnsfp_vep.py"
_SPEC = importlib.util.spec_from_file_location("prepare_dbnsfp_vep", _SCRIPT)
assert _SPEC and _SPEC.loader
_MOD = importlib.util.module_from_spec(_SPEC)
sys.modules["prepare_dbnsfp_vep"] = _MOD
_SPEC.loader.exec_module(_MOD)

_chrom_from_path = _MOD._chrom_from_path
_load_panel_intervals = _MOD._load_panel_intervals
_main = _MOD.main


def test_chrom_from_path_parses_academic_names():
    assert _chrom_from_path(Path("dbNSFP5.4a_variant.chr1.gz")) == "1"
    assert _chrom_from_path(Path("dbNSFP5.4a_variant.chr22.gz")) == "22"
    assert _chrom_from_path(Path("dbNSFP5.4a_variant.chrX.gz")) == "X"
    assert _chrom_from_path(Path("dbNSFP5.4a_variant.chrY.gz")) == "Y"
    assert _chrom_from_path(Path("dbNSFP5.4a_variant.chrM.gz")) == "M"


def test_panel_intervals_strip_chr_prefix(tmp_path: Path):
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t99\t200\tGENE\nchrX\t0\t10\tXGENE\n", encoding="utf-8")
    intervals = _load_panel_intervals(bed)
    assert intervals is not None
    assert intervals["1"] == [(100, 200)]
    assert intervals["X"] == [(1, 10)]


def test_prepare_subsets_to_panel_and_writes_tabix(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    if not shutil.which("bgzip") or not shutil.which("tabix"):
        pytest.skip("bgzip/tabix not installed")

    header = "#chr\tpos(1-based)\tref\talt\tCADD_phred\tREVEL_score\tfathmm-XF_coding_pred\n"
    rows = [
        "1\t50\tA\tG\t1.0\t0.1\tN\n",
        "1\t150\tC\tT\t28.5\t0.85\tD\n",
        "1\t250\tG\tA\t2.0\t0.2\tN\n",
    ]
    chrom = tmp_path / "dbNSFP5.4a_variant.chr1.gz"
    with gzip.open(chrom, "wt", encoding="utf-8") as handle:
        handle.write(header)
        handle.writelines(rows)
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t100\t200\tPINK1\n", encoding="utf-8")
    out_prefix = tmp_path / "dbNSFP5.4a_panel"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prepare_dbnsfp_vep.py",
            "--out-prefix",
            str(out_prefix),
            "--panel-bed",
            str(bed),
            str(chrom),
        ],
    )
    assert _main() == 0
    out_gz = Path(str(out_prefix) + ".gz")
    assert out_gz.is_file()
    assert Path(str(out_gz) + ".tbi").is_file()
    with gzip.open(out_gz, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert lines[0].startswith("#chr")
    assert len(lines) == 2
    assert "28.5" in lines[1]
    assert "0.85" in lines[1]
