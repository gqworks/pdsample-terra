"""Lift gnomAD AF from dbNSFP CSQ into INFO."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "lift_gnomad_from_csq.py"
_SPEC = importlib.util.spec_from_file_location("lift_gnomad_from_csq", _SCRIPT)
assert _SPEC and _SPEC.loader
_MOD = importlib.util.module_from_spec(_SPEC)
sys.modules["lift_gnomad_from_csq"] = _MOD
_SPEC.loader.exec_module(_MOD)

lift = _MOD.lift


def test_lift_copies_joint_af_from_csq(tmp_path: Path):
    csq = "Allele|Consequence|gnomAD4.1_joint_AF|gnomAD4.1_joint_POPMAX_AF"
    src = tmp_path / "in.vcf"
    src.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {csq}">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                f"1\t100\t.\tA\tG\t.\tPASS\tCSQ=G|missense_variant|0.001|0.002",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    lift(src, out)
    text = out.read_text(encoding="utf-8")
    assert "##INFO=<ID=gnomAD_AF," in text
    assert "##INFO=<ID=AF_popmax," in text
    assert "gnomAD_AF=0.001" in text
    assert "AF_popmax=0.002" in text
