"""Annotate sample VCF with panel-sliced ClinVar P/LP sites."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "annotate_clinvar_plp.py"
_SPEC = importlib.util.spec_from_file_location("annotate_clinvar_plp", _SCRIPT)
assert _SPEC and _SPEC.loader
_MOD = importlib.util.module_from_spec(_SPEC)
sys.modules["annotate_clinvar_plp"] = _MOD
_SPEC.loader.exec_module(_MOD)

annotate = _MOD.annotate


def test_annotate_plp_allele_and_sv_id(tmp_path: Path):
    sites = tmp_path / "sites.vcf"
    sites.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr6\t161350208\t.\tC\tT\t.\tPASS\tCLNSIG=Pathogenic",
                "chr6\t161569000\tdelPRKN\tN\t<DEL>\t.\tPASS\tCLNSIG=Pathogenic",
                "chr1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Benign",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample = tmp_path / "sample.vcf"
    sample.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr6\t161350208\t.\tC\tT\t.\tPASS\t.",
                "chr6\t161569000\tdelPRKN\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL",
                "chr1\t100\t.\tA\tG\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    matched = annotate(sample, sites, out)
    text = out.read_text(encoding="utf-8")
    assert matched == 2
    assert "##INFO=<ID=ClinVar_CLNSIG," in text
    assert "ClinVar_CLNSIG=Pathogenic" in text
    assert "chr1\t100\t.\tA\tG\t.\tPASS\t." in text
