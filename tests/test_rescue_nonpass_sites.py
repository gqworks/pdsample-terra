"""Restore non-PASS ClinVar P/LP and GBA1 browser sites."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "rescue_nonpass_sites.py"
_SPEC = importlib.util.spec_from_file_location("rescue_nonpass_sites", _SCRIPT)
assert _SPEC and _SPEC.loader
_MOD = importlib.util.module_from_spec(_SPEC)
sys.modules["rescue_nonpass_sites"] = _MOD
_SPEC.loader.exec_module(_MOD)

rescue = _MOD.rescue


def test_rescue_clinvar_and_gba1_nonpass(tmp_path: Path):
    clinvar = tmp_path / "clinvar.vcf"
    clinvar.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr6\t161350208\t.\tC\tT\t.\tPASS\tCLNSIG=Pathogenic",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    browser = tmp_path / "gba1.tsv"
    browser.write_text(
        "cdna\tprotein\tvariant_name\trsid\tpd_class\tgd_acmg\tpos_hg38\tallele_frequency\treferences\tcomments\n"
        "c.635C>G\tp.(Ser212Ter)\tp.S212*\t.\tsevere\tPathogenic\t155238246\t0\t1\t\n",
        encoding="utf-8",
    )
    sample = tmp_path / "sample.vcf"
    sample.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                "chr6\t161350208\t.\tC\tT\t.\thard_filter\t.\tGT\t0/1",
                "chr1\t155238246\t.\tC\tG\t.\thard_filter\t.\tGT\t0/1",
                "chr2\t100\t.\tA\tT\t.\thard_filter\t.\tGT\t0/1",
                "chr6\t161350208\t.\tC\tA\t.\thard_filter\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    count = rescue(sample, out, clinvar_vcf=clinvar, gba1_browser=browser)
    text = out.read_text(encoding="utf-8")
    assert count == 2
    assert "PDVAR_RESCUE" in text
    assert "chr6\t161350208\t.\tC\tT\t.\tPASS\tPDVAR_RESCUE\tGT\t0/1" in text
    assert "chr1\t155238246\t.\tC\tG\t.\tPASS\tPDVAR_RESCUE\tGT\t0/1" in text
    assert "chr2\t100\t.\tA\tT\t.\thard_filter\t.\tGT\t0/1" in text
    assert "chr6\t161350208\t.\tC\tA\t.\thard_filter\t.\tGT\t0/1" in text
