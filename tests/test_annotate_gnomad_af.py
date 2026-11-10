"""Annotate sample VCF with panel-sliced gnomAD joint AF."""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "annotate_gnomad_af.py"
_SPEC = importlib.util.spec_from_file_location("annotate_gnomad_af", _SCRIPT)
assert _SPEC and _SPEC.loader
_MOD = importlib.util.module_from_spec(_SPEC)
sys.modules["annotate_gnomad_af"] = _MOD
_SPEC.loader.exec_module(_MOD)

annotate = _MOD.annotate


def test_annotate_intronic_joint_af(tmp_path: Path):
    """ADAR-style intron allele gets AF_joint even though dbNSFP would miss it."""
    sites = tmp_path / "sites.vcf"
    sites.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=AF_joint,Number=A,Type=Float,Description="joint AF">',
                '##INFO=<ID=AF_grpmax_joint,Number=A,Type=Float,Description="grpmax">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t154588785\t.\tG\tC\t.\tPASS\tAF_joint=0.9852;AF_grpmax_joint=0.9997",
                "chr1\t226894490\t.\tG\tA\t.\tPASS\tAF_joint=0.00645;AF_grpmax_joint=0.01",
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
                "chr1\t154588785\t.\tG\tC\t.\tPASS\tCSQ=C|intron_variant",
                "chr1\t226894490\t.\tG\tA\t.\tPASS\t.",
                "chr1\t100\t.\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    matched = annotate(sample, sites, out)
    text = out.read_text(encoding="utf-8")
    assert matched == 2
    assert "##INFO=<ID=gnomAD_AF," in text
    assert "##INFO=<ID=AF_popmax," in text
    assert "gnomAD_AF=0.9852" in text
    assert "AF_popmax=0.9997" in text
    assert "gnomAD_AF=0.00645" in text
    assert "chr1\t100\t.\tA\tT\t.\tPASS\t." in text


def test_annotate_accepts_bare_chrom_and_standard_af(tmp_path: Path):
    sites = tmp_path / "sites.vcf"
    sites.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
                '##INFO=<ID=AF_grpmax,Number=A,Type=Float,Description="grpmax">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "1\t100\t.\tA\tG\t.\tPASS\tAF=0.02;AF_grpmax=0.03",
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
                "chr1\t100\t.\tA\tG\t.\tPASS\tDP=10",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    annotate(sample, sites, out)
    text = out.read_text(encoding="utf-8")
    assert "gnomAD_AF=0.02" in text
    assert "AF_popmax=0.03" in text
    assert "DP=10" in text


def test_annotate_reads_bgz_sites(tmp_path: Path):
    """Panel sites are staged as .vcf.bgz; must not be opened as plain UTF-8."""
    sites_text = (
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=AF_joint,Number=A,Type=Float,Description="joint AF">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t154588785\t.\tG\tC\t.\tPASS\tAF_joint=0.9852",
            ]
        )
        + "\n"
    )
    sites = tmp_path / "sites.vcf.bgz"
    with gzip.open(sites, "wt", encoding="utf-8") as handle:
        handle.write(sites_text)
    sample = tmp_path / "sample.vcf.gz"
    with gzip.open(sample, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                    "chr1\t154588785\t.\tG\tC\t.\tPASS\t.",
                ]
            )
            + "\n"
        )
    out = tmp_path / "out.vcf"
    matched = annotate(sample, sites, out)
    text = out.read_text(encoding="utf-8")
    assert matched == 1
    assert "gnomAD_AF=0.9852" in text
