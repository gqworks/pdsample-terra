"""Empty / missing / zero-byte Terra priors must not count as reusable files."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from importlib.machinery import SourceFileLoader

ROOT = Path(__file__).resolve().parents[1]
probe = SourceFileLoader(
    "probe_prior_file",
    str(ROOT / "workflow" / "scripts" / "probe_prior_file.py"),
).load_module()


def test_blank_and_placeholder_paths_are_unset():
    for path in ("", "  ", "null", "None", "NA", ".", "undefined"):
        assert probe.usable(path, kind="introme") is False
        assert probe.usable(path, kind="vcf") is False


def test_missing_and_zero_byte_files_are_unset(tmp_path: Path):
    missing = tmp_path / "gone.tsv"
    empty = tmp_path / "empty.tsv"
    empty.write_bytes(b"")
    assert probe.usable(str(missing), kind="introme") is False
    assert probe.usable(str(empty), kind="introme") is False
    assert probe.usable(str(empty), kind="vcf") is False
    assert probe.usable(str(empty), kind="any") is False


def test_real_introme_tsv_is_usable(tmp_path: Path):
    path = tmp_path / "APGS.introme.predictions.tsv"
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tgene\tsamples_with_alt\tintrome_score\nchr1\t1\tA\tT\tPRKN\t1\t0.91\n",
        encoding="utf-8",
    )
    assert probe.usable(str(path), kind="introme") is True
    assert probe.usable(str(path), kind="vcf") is False


def test_placeholder_text_is_not_usable(tmp_path: Path):
    path = tmp_path / "this.introme"
    path.write_text("placeholder\n", encoding="utf-8")
    assert probe.usable(str(path), kind="introme") is False
    assert probe.usable(str(path), kind="any") is False
    assert probe.usable(str(path), kind="json") is False
    assert probe.usable(str(path), kind="tsv") is False


def test_json_csv_tsv_kinds(tmp_path: Path):
    json_path = tmp_path / "sample.gauchian.json"
    json_path.write_text('{"Sample": "APGS", "CopyNumber": 2}\n', encoding="utf-8")
    csv_path = tmp_path / "panel.csv"
    csv_path.write_text("Gene_ID,Symbol,inheritance_primary\nENSG1,GBA1,AD\n", encoding="utf-8")
    tsv_path = tmp_path / "gba1.tsv"
    tsv_path.write_text("cdna\tprotein\tpos_hg38\nc.1A>G\tp.M1V\t1\n", encoding="utf-8")
    avi_path = tmp_path / "avi.tsv"
    avi_path.write_text("chrom\tpos\tref\talt\tavi_raw\tavi_quantile\tavi_phred\tavi_url\n", encoding="utf-8")
    assert probe.usable(str(json_path), kind="gauchian_json") is True
    assert probe.usable(str(json_path), kind="vcf") is False
    assert probe.usable(str(csv_path), kind="panel_csv") is True
    assert probe.usable(str(csv_path), kind="tsv") is False
    assert probe.usable(str(tsv_path), kind="gba1_pd_browser") is True
    assert probe.usable(str(avi_path), kind="avi_tsv") is True
    assert probe.usable(str(avi_path), kind="vcf") is False


def test_real_vcf_is_usable(tmp_path: Path):
    vcf = tmp_path / "sample.vep.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )
    gz = tmp_path / "sample.vep.vcf.gz"
    gz.write_bytes(gzip.compress(vcf.read_bytes()))
    assert probe.usable(str(vcf), kind="vcf") is True
    assert probe.usable(str(gz), kind="vep_vcf") is True
    assert probe.usable(str(vcf), kind="introme") is False


def test_cli_writes_usable_txt(tmp_path: Path):
    vcf = tmp_path / "x.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    out = tmp_path / "usable.txt"
    assert probe.main(["--path", str(vcf), "--kind", "vcf", "--out", str(out)]) == 0
    assert out.read_text(encoding="utf-8").strip() == "true"
    empty = tmp_path / "z"
    empty.write_bytes(b"")
    assert probe.main(["--path", str(empty), "--kind", "vcf", "--out", str(out)]) == 0
    assert out.read_text(encoding="utf-8").strip() == "false"
    assert probe.main(["--path", "", "--kind", "introme", "--out", str(out)]) == 0
    assert out.read_text(encoding="utf-8").strip() == "false"
