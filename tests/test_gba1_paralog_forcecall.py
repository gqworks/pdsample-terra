"""Paralog-aware GBA1 hotspot force-call helpers and VCF ingest."""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

from pdvar import callsets

ROOT = Path(__file__).resolve().parents[1]
_SCRIPT = ROOT / "workflow" / "scripts" / "gba1_paralog_forcecall.py"
_SPEC = importlib.util.spec_from_file_location("gba1_paralog_forcecall", _SCRIPT)
assert _SPEC and _SPEC.loader
_MOD = importlib.util.module_from_spec(_SPEC)
sys.modules["gba1_paralog_forcecall"] = _MOD
_SPEC.loader.exec_module(_MOD)

DISC = ROOT / "resources" / "gba1" / "gba1_homology_discriminators.tsv"
HOT = ROOT / "resources" / "gba1" / "gba1_homology_hotspots.tsv"
PATCH = ROOT / "resources" / "gba1" / "gbap1_grch38_patches.tsv"


def test_load_shipped_assets():
    discs = _MOD.load_discriminators(DISC)
    hotspots = _MOD.load_hotspots(HOT)
    patches = _MOD.load_patches(PATCH)
    assert any(d.label == "L483P" for d in discs)
    assert {h.token for h in hotspots} >= {"L483P", "A495P", "D448H"}
    assert {p.pos for p in patches} >= {155214266, 155214276, 155214576, 155214590}
    l483 = next(h for h in hotspots if h.token == "L483P")
    assert l483.cdna == "c.1448T>C"
    assert l483.pos == 155235252
    assert l483.ref == "A" and l483.alt == "G"


def test_apply_patches_restores_major_gbap1_alleles():
    start = 155214266
    length = 155214590 - start + 1
    chars = ["A"] * length
    chars[155214266 - start] = "A"
    chars[155214276 - start] = "C"
    chars[155214576 - start] = "C"
    chars[155214590 - start] = "C"
    raw = "".join(chars)
    patches = _MOD.load_patches(PATCH)
    patched, applied = _MOD.apply_patches(raw, window_start=start, patches=patches, chrom="chr1")
    assert len(applied) == 4
    assert patched[155214266 - start] == "G"
    assert patched[155214276 - start] == "T"
    assert patched[155214576 - start] == "G"
    assert patched[155214590 - start] == "G"


def test_template_score_prefers_gene_for_l483p_after_patch():
    """Orimo Case 1: L483P read matches unpatched GBAP1; prefers gene after patch."""
    # Longer asymmetric flanks so reverse-complement does not tie scores.
    flank_l, flank_r = "ATGCCGTAA", "TTGACATCG"
    gene = flank_l + "C" + "C" + "A" + flank_r
    gene_l483p = flank_l + "C" + "C" + "G" + flank_r
    grch38_pseudo = flank_l + "C" + "C" + "G" + flank_r
    patched_pseudo = flank_l + "G" + "G" + "G" + flank_r

    g_unpatched, p_unpatched = _MOD.score_read_vs_templates(gene_l483p, gene, grch38_pseudo)
    assert p_unpatched > g_unpatched

    g_patched, p_patched = _MOD.score_read_vs_templates(gene_l483p, gene, patched_pseudo)
    assert g_patched > p_patched


def test_genotype_l444p_from_gene_assigned_counts():
    gt, dp, ad_ref, ad_alt = _MOD.genotype_from_counts(
        Counter({"A": 20, "G": 8}),
        ref="A",
        alt="G",
        min_alt_reads=3,
        min_alt_fraction=0.15,
        min_gene_dp=8,
    )
    assert gt == "0/1"
    assert dp == 28 and ad_ref == 20 and ad_alt == 8
    gt0, *_ = _MOD.genotype_from_counts(
        Counter({"A": 30}),
        ref="A",
        alt="G",
        min_alt_reads=3,
        min_alt_fraction=0.15,
        min_gene_dp=8,
    )
    assert gt0 == "0/0"


def test_write_vcf_and_parse_gba1_paralog(tmp_path: Path):
    sites = [
        {
            "chrom": "chr1",
            "pos": 155235252,
            "ref": "A",
            "alt": "G",
            "token": "L483P",
            "cdna": "c.1448T>C",
            "protein": "p.Leu483Pro",
            "pd_class": "severe",
            "gt": "0/1",
            "dp": 28,
            "ad_ref": 20,
            "ad_alt": 8,
        },
        {
            "chrom": "chr1",
            "pos": 155235217,
            "ref": "C",
            "alt": "G",
            "token": "A495P",
            "cdna": "c.1483G>C",
            "protein": "p.Ala495Pro",
            "pd_class": "severe",
            "gt": "0/0",
            "dp": 40,
            "ad_ref": 40,
            "ad_alt": 0,
        },
    ]
    vcf = tmp_path / "sample.gba1_paralog.vcf"
    _MOD.write_vcf(vcf, "MONOPDAUS_000046", sites)
    text = vcf.read_text(encoding="utf-8")
    assert "PDVAR_GBA1_PARALOG" in text
    assert "TOKEN=L483P" in text
    table, seen, carriers = callsets.parse_gba1_paralog(vcf)
    assert seen == {"MONOPDAUS_000046"}
    assert carriers == {"MONOPDAUS_000046"}
    assert len(table) == 1
    row = table.iloc[0]
    assert row["hit_source"] == "gba1_paralog"
    assert row["hgvsc"] == "c.1448T>C"
    assert row["gene"] == "GBA1"
    assert row["clinvar_class"] == "P/LP"
    assert row["depth"] == 28.0


def test_drop_redundant_prefers_gauchian_over_paralog():
    df = pd.DataFrame(
        [
            {
                "gene": "GBA1",
                "hit_source": "gauchian",
                "hgvsc": "c.1448T>C",
                "dosage_event": "other_variant",
                "variant_id": "g",
                "gba1_pd_class": "severe",
            },
            {
                "gene": "GBA1",
                "hit_source": "gba1_paralog",
                "hgvsc": "c.1448T>C",
                "dosage_event": "paralog_hotspot",
                "variant_id": "p",
                "gba1_pd_class": "severe",
            },
        ]
    )
    out = callsets.drop_redundant_gauchian_snvs(df)
    assert list(out["hit_source"]) == ["gauchian"]


@pytest.mark.skipif(importlib.util.find_spec("pysam") is None, reason="pysam not installed")
def test_forcecall_patched_template_recovers_l444p(tmp_path: Path):
    """Read on GBAP1 coords with gene-like sequence; patched templates assign to gene."""
    import pysam

    disc_path = tmp_path / "disc.tsv"
    disc_path.write_text(
        "chrom\tgene_pos\tgene_allele\tpseudo_pos\tpseudo_allele\tlabel\n"
        "chr1\t1000\tC\t2000\tT\tanchor\n"
        "chr1\t1100\tC\t2100\tG\tV499=\n"
        "chr1\t1110\tC\t2110\tG\tA495P\n"
        "chr1\t1120\tA\t2120\tG\tL483P\n"
        "chr1\t1200\tA\t2200\tG\tanchor\n",
        encoding="utf-8",
    )
    hot_path = tmp_path / "hot.tsv"
    hot_path.write_text(
        "token\tcdna\tprotein\tchrom\tpos\tref\talt\tpseudo_pos\tpd_class\n"
        "L483P\tc.1448T>C\tp.Leu483Pro\tchr1\t1120\tA\tG\t2120\tsevere\n",
        encoding="utf-8",
    )
    patch_path = tmp_path / "patch.tsv"
    patch_path.write_text(
        "chrom\tpos\tgrch38_allele\tpatched_allele\tnote\n"
        "chr1\t2000\tC\tT\tanchor\n"
        "chr1\t2100\tC\tG\tV499=\n"
        "chr1\t2110\tC\tG\tA495P\n",
        encoding="utf-8",
    )

    chrom = ["N"] * 2300

    def put(pos: int, base: str) -> None:
        chrom[pos - 1] = base

    for pos, base in [(1000, "C"), (1100, "C"), (1110, "C"), (1120, "A"), (1200, "A")]:
        put(pos, base)
    for i in range(1000, 1201):
        if chrom[i - 1] == "N":
            chrom[i - 1] = "A"
    for pos, base in [(2000, "C"), (2100, "C"), (2110, "C"), (2120, "G"), (2200, "G")]:
        put(pos, base)
    for i in range(2000, 2201):
        if chrom[i - 1] == "N":
            chrom[i - 1] = "T"

    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "".join(chrom) + "\n", encoding="utf-8")
    pysam.faidx(str(fa))

    gene_slice = list("".join(chrom[999:1200]))
    gene_slice[1120 - 1000] = "G"
    read_seq = "".join(gene_slice)
    bam_path = tmp_path / "misaligned.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": len(chrom)}],
        "RG": [{"ID": "rg1", "SM": "MH06245_SYD"}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        read = pysam.AlignedSegment()
        read.query_name = "gene_frag_on_pseudo"
        read.query_sequence = read_seq
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 1999
        read.mapping_quality = 0
        read.cigar = [(0, len(read_seq))]
        read.next_reference_id = -1
        read.next_reference_start = -1
        read.template_length = 0
        read.query_qualities = pysam.qualitystring_to_array("I" * len(read_seq))
        out.write(read)
    pysam.index(str(bam_path))

    discs = _MOD.load_discriminators(disc_path)
    hotspots = _MOD.load_hotspots(hot_path)
    patches = _MOD.load_patches(patch_path)
    sites, summary = _MOD.forcecall(
        bam_path,
        reference=fa,
        discriminators=discs,
        hotspots=hotspots,
        patches=patches,
        sample_id="MH06245_SYD",
        min_alt_reads=1,
        min_alt_fraction=0.1,
        min_gene_dp=1,
        mapq_floor=0,
        pad=0,
    )
    assert summary["assignment_mode"] == "patched_grch38_template"
    assert summary["gene_assigned_fragments"] >= 1
    assert summary["templates"]["pseudo_patched"]["patches_applied"]
    l483 = next(s for s in sites if s["token"] == "L483P")
    assert l483["gt"] in {"0/1", "1/1"}
    assert l483["ad_alt"] >= 1
