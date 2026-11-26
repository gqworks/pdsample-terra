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


def test_load_shipped_assets():
    discs = _MOD.load_discriminators(DISC)
    hotspots = _MOD.load_hotspots(HOT)
    assert any(d.label == "L483P" for d in discs)
    assert {h.token for h in hotspots} >= {"L483P", "A495P", "D448H"}
    l483 = next(h for h in hotspots if h.token == "L483P")
    assert l483.cdna == "c.1448T>C"
    assert l483.pos == 155235252
    assert l483.ref == "A" and l483.alt == "G"


def test_assign_locus_from_discriminator_bases():
    discs = _MOD.load_discriminators(DISC)
    hotspots = _MOD.load_hotspots(HOT)
    scoring = _MOD.discriminators_for_scoring(discs, hotspots)
    assert {d.label for d in scoring} == {"anchor", "V499="}
    assert len(scoring) == 5
    # Gene-like bases at gene coordinates (anchors + V499=).
    gene_obs = {
        ("1", 155234903): "C",
        ("1", 155235203): "C",
        ("1", 155235379): "A",
        ("1", 155235412): "G",
        ("1", 155235918): "A",
    }
    g, p = _MOD.score_assignment(gene_obs, scoring, locus="gene")
    assert g > p
    assert _MOD.assign_locus(g, p) == "gene"
    # Pseudogene-like bases seen on a GBAP1-mapped fragment (anchors + V499=).
    pseudo_obs = {
        ("1", 155214276): "T",
        ("1", 155214576): "G",
        ("1", 155214752): "G",
        ("1", 155214785): "A",
        ("1", 155215236): "G",
    }
    g2, p2 = _MOD.score_assignment(pseudo_obs, scoring, locus="pseudo")
    assert p2 > g2
    assert _MOD.assign_locus(g2, p2) == "pseudo"


def test_hotspot_alts_do_not_vote_fragments_toward_pseudo():
    """RecNciI-like fragment: gene anchors + GBAP1 alleles at L483P/A495P/D448H.

    Before excluding hotspots from scoring, those pathogenic alts tip the
    fragment to ``pseudo`` and empty the gene pileup (MH06245_SYD pattern).
    """
    discs = _MOD.load_discriminators(DISC)
    hotspots = _MOD.load_hotspots(HOT)
    scoring = _MOD.discriminators_for_scoring(discs, hotspots)
    obs = {
        ("1", 155234903): "C",  # gene anchor
        ("1", 155235203): "C",  # V499= gene
        ("1", 155235217): "G",  # A495P pathogenic (= pseudo allele)
        ("1", 155235252): "G",  # L483P pathogenic (= pseudo allele)
        ("1", 155235727): "G",  # D448H pathogenic (= pseudo allele)
        ("1", 155235379): "A",  # gene anchor
    }
    g_all, p_all = _MOD.score_assignment(obs, discs, locus="gene")
    # Hotspot alts cancel gene anchors → not assigned to gene (tie or pseudo).
    assert _MOD.assign_locus(g_all, p_all) != "gene"
    g, p = _MOD.score_assignment(obs, scoring, locus="gene")
    assert g > p
    assert _MOD.assign_locus(g, p) == "gene"


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
def test_forcecall_recovers_l444p_from_gbap1_mapped_read(tmp_path: Path):
    """Read mapped to GBAP1 but carrying gene-like flanking SNPs + L483P alt."""
    import pysam

    bam_path = tmp_path / "misaligned.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 249250621}],
        "RG": [{"ID": "rg1", "SM": "MONOPDAUS_000046"}],
    }
    # Cover GBAP1 images of exon 9–11 discriminators. Gene-like bases at those
    # coords assign the fragment to GBA1; G at 155214625 is pathogenic L483P.
    start = 155214270  # 0-based → covers 155214276 .. 155215236
    end_1based = 155215250
    length = end_1based - start
    seq = ["N"] * length

    def set_base(pos_1based: int, base: str) -> None:
        idx = pos_1based - (start + 1)
        assert 0 <= idx < length, pos_1based
        seq[idx] = base

    set_base(155214276, "C")
    set_base(155214576, "C")
    set_base(155214590, "C")
    set_base(155214625, "G")  # L483P alt on gene-assigned fragment
    set_base(155214752, "A")
    set_base(155214785, "G")
    set_base(155215100, "C")
    set_base(155215236, "A")
    read_seq = "".join(seq)
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        read = pysam.AlignedSegment()
        read.query_name = "gene_frag_on_pseudo"
        read.query_sequence = read_seq
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 0
        read.cigar = [(0, length)]
        read.next_reference_id = -1
        read.next_reference_start = -1
        read.template_length = 0
        read.query_qualities = pysam.qualitystring_to_array("I" * length)
        out.write(read)
    pysam.index(str(bam_path))

    discs = _MOD.load_discriminators(DISC)
    hotspots = _MOD.load_hotspots(HOT)
    sites, summary = _MOD.forcecall(
        bam_path,
        discriminators=discs,
        hotspots=hotspots,
        sample_id="MONOPDAUS_000046",
        min_alt_reads=1,
        min_alt_fraction=0.1,
        min_gene_dp=1,
        mapq_floor=0,
    )
    assert summary["gene_assigned_fragments"] >= 1
    l483 = next(s for s in sites if s["token"] == "L483P")
    assert l483["gt"] in {"0/1", "1/1"}
    assert l483["ad_alt"] >= 1
