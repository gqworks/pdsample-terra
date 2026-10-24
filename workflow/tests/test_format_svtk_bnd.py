"""Breakends must be dictionary-ordered before GATK SVCluster."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflow" / "wdl" / "vendor" / "gatk-sv" / "format_svtk_vcf_for_gatk.py"


def _load():
    spec = importlib.util.spec_from_file_location("format_svtk_vcf_for_gatk", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_breakend_with_mate_before_pos_is_flipped():
    fmt = _load()
    pos, end2, strands, alt = fmt.order_same_chrom_breakend(5000, 1000, "+-", "N", "chr1")
    assert (pos, end2, strands) == (1000, 5000, "-+")
    assert alt == "]chr1:5000]N"


def test_already_ordered_breakend_is_unchanged():
    fmt = _load()
    pos, end2, strands, alt = fmt.order_same_chrom_breakend(1000, 5000, "+-", "N", "chr1")
    assert (pos, end2, strands, alt) == (1000, 5000, "+-", None)


def test_breakend_alt_matches_gatk_strand_encoding():
    fmt = _load()
    assert fmt.breakend_alt("N", "chr2", 40, "+-") == "N[chr2:40["
    assert fmt.breakend_alt("N", "chr2", 40, "++") == "N]chr2:40]"
    assert fmt.breakend_alt("N", "chr2", 40, "-+") == "]chr2:40]N"
    assert fmt.breakend_alt("N", "chr2", 40, "--") == "[chr2:40[N"


def test_interchrom_breakend_moves_to_earlier_contig():
    fmt = _load()
    order = {"chr1": 0, "chr2": 1, "chr3": 2}
    # Record on chr2 with mate on chr1 must relocate to chr1.
    chrom, pos, chr2, end2, strands, alt = fmt.order_breakends(
        "chr2", 5000, "chr1", 1000, "+-", "N", order
    )
    assert (chrom, pos, chr2, end2, strands) == ("chr1", 1000, "chr2", 5000, "-+")
    assert alt == "]chr2:5000]N"


def test_interchrom_already_ordered_is_unchanged():
    fmt = _load()
    order = {"chr1": 0, "chr2": 1}
    chrom, pos, chr2, end2, strands, alt = fmt.order_breakends(
        "chr1", 1000, "chr2", 5000, "+-", "N", order
    )
    assert (chrom, pos, chr2, end2, strands, alt) == ("chr1", 1000, "chr2", 5000, "+-", None)


def test_contig_order_uses_fai_not_lexicographic():
    """chr10 must not sort before chr2 (lexicographic trap)."""
    fmt = _load()
    order = {"chr1": 0, "chr2": 1, "chr10": 9}
    # chr10:100 with mate chr2:5000 — chr2 precedes chr10 in dictionary order.
    chrom, pos, chr2, end2, strands, alt = fmt.order_breakends(
        "chr10", 100, "chr2", 5000, "++", "N", order
    )
    assert (chrom, pos, chr2, end2) == ("chr2", 5000, "chr10", 100)
    assert strands == "++"
    assert alt == "N]chr10:100]"


def test_load_contig_order_from_genome_file(tmp_path):
    fmt = _load()
    genome = tmp_path / "genome.file"
    genome.write_text("chr1\t248956422\nchr2\t242193529\nchr10\t133797422\n", encoding="utf-8")
    assert fmt.load_contig_order(genome) == {"chr1": 0, "chr2": 1, "chr10": 2}
