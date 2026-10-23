"""Same-chromosome BNDs must be ordered before GATK SVCluster."""

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
