"""Figure 3 and Figure 4 splice-heuristic calls on a tiny reference."""

import pandas as pd
import pytest

from pdvar.splice_heuristics import (
    Exon,
    IndexedFasta,
    MemoryFasta,
    annotate_splice_heuristics,
    load_gtf_exons,
)

# Plus-strand exon, 1-based bases 201-400. Donor +1 is 401. Acceptor E+1 is 201.
PLUS = Exon("1", 200, 400, "+", "GENE", "ENST0001")
MINUS = Exon("1", 600, 700, "-", "GENE", "ENST0009")


def _call(pos, ref, alt, edits=None, exons=None):
    seq = ["A"] * 800
    for genomic, base in (edits or {}).items():
        seq[genomic - 1] = base
    seq[pos - 1] = ref
    fasta = MemoryFasta({"1": "".join(seq)})
    frame = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": "Gene",
                "hgvsc": "ENST0001.1:c.10A>G",
                "evidence_score": 4,
            }
        ]
    )
    return annotate_splice_heuristics(frame, exons if exons is not None else [PLUS], fasta).iloc[0]


def _expect(row, rule, subgroup, percent, count):
    assert row["splice_rule"] == rule
    assert row["splice_subgroup"] == subgroup
    assert row["splice_sav_pct"] == pytest.approx(percent)
    assert int(row["splice_heuristic_n"]) == count
    assert row["evidence_score"] == 4


# pos, ref, alt, fasta edits, rule, subgroup, %SAV, n
CASES = [
    ("DD2 standard", 401, "G", "A", {402: "T"}, "DD2", "standard", 99.9, 2218),
    ("DD1 context", 402, "C", "T", {401: "G"}, "DD1", "context", 81.2, 16),
    (
        "DD1 context ideal",
        402,
        "C",
        "A",
        {398: "C", 400: "G", 401: "G", 405: "G", 406: "T"},
        "DD1",
        "context",
        81.2,
        16,
    ),
    ("DD1 standard", 401, "G", "A", {402: "C"}, "DD1", "standard", 100.0, 28),
    ("DD3 standard", 398, "C", "T", {}, "DD3", "standard", 100.0, 78),
    ("DD12 standard", 398, "C", "T", {399: "G"}, "DD12", "standard", 68.2, 22),
    ("DD12 auxiliary", 398, "A", "C", {}, "DD12", "auxiliary", 45.3, 64),
    ("DD11 standard", 399, "A", "G", {}, "DD11", "standard", 91.1, 146),
    ("DD11 auxiliary", 399, "G", "C", {}, "DD11", "auxiliary", 38.6, 57),
    ("DD4 standard", 405, "G", "A", {}, "DD4", "standard", 99.7, 583),
    ("DD7 standard", 400, "G", "A", {405: "G"}, "DD7", "standard", 94.0, 667),
    ("DD7 auxiliary", 400, "A", "G", {405: "G"}, "DD7", "auxiliary", 25.0, 12),
    ("DD6 auxiliary", 405, "A", "C", {}, "DD6", "auxiliary", 28.6, 14),
    ("DD6 standard", 405, "G", "C", {400: "G"}, "DD6", "standard", 95.1, 446),
    ("DD5 standard", 403, "A", "C", {404: "C"}, "DD5", "standard", 95.4, 108),
    ("DD8 standard", 403, "A", "C", {}, "DD8", "standard", 93.6, 344),
    ("DD8 auxiliary", 403, "C", "A", {}, "DD8", "auxiliary", 25.0, 20),
    ("DD10 standard", 404, "A", "G", {}, "DD10", "standard", 91.5, 165),
    ("DD10 auxiliary", 404, "G", "C", {}, "DD10", "auxiliary", 14.3, 28),
    ("DD9 standard", 406, "T", "A", {}, "DD9", "standard", 91.7, 120),
    ("DD9 auxiliary", 406, "A", "T", {}, "DD9", "auxiliary", 19.0, 63),
    ("DA2 standard", 200, "G", "A", {199: "A"}, "DA2", "standard", 99.7, 2817),
    ("DA1 standard", 191, "C", "A", {192: "G"}, "DA1", "standard", 97.7, 304),
    ("DA1 context", 181, "C", "A", {182: "G", 159: "T", 161: "A"}, "DA1", "context", 78.0, 100),
    ("DA1 auxiliary", 196, "C", "A", {197: "G"}, "DA1", "auxiliary", 28.0, 232),
    ("DA6 standard", 201, "G", "A", {198: "C"}, "DA6", "standard", 66.7, 96),
    ("DA6 auxiliary", 201, "A", "G", {198: "C"}, "DA6", "auxiliary", 14.8, 27),
    ("DA3 context", 201, "G", "A", {198: "T"}, "DA3", "context", 70.4, 81),
    ("DA3 standard", 201, "G", "A", {}, "DA3", "standard", 99.2, 124),
    ("DA4 context", 198, "C", "T", {201: "G"}, "DA4", "context", 23.8, 21),
    ("DA4 standard", 198, "C", "A", {201: "G"}, "DA4", "standard", 96.8, 157),
    ("DA4 auxiliary", 198, "T", "C", {}, "DA4", "auxiliary", 17.4, 23),
    ("DA8 context", 197, "T", "A", {genomic: "T" for genomic in range(177, 198)}, "DA8", "context", 7.7, 39),
    ("DA8 standard", 197, "A", "C", {}, "DA8", "standard", 33.3, 57),
    ("DA5 standard", 196, "T", "A", {}, "DA5", "standard", 93.1, 102),
    ("DA5 auxiliary", 196, "T", "C", {}, "DA5", "auxiliary", 43.3, 164),
    ("DA9 context", 191, "T", "A", {genomic: "T" for genomic in range(177, 198)}, "DA9", "context", 3.3, 153),
    ("DA9 auxiliary", 191, "A", "T", {}, "DA9", "auxiliary", 4.6, 524),
    ("DA9 standard", 191, "C", "A", {}, "DA9", "standard", 25.4, 508),
    ("DA7 standard", 159, "T", "A", {161: "A"}, "DA7", "standard", 38.9, 18),
    ("DA7 context", 159, "T", "A", {161: "A", 169: "T", 171: "A"}, "DA7", "context", 19.2, 151),
    ("DA7 auxiliary", 160, "A", "C", {159: "T", 161: "A"}, "DA7", "auxiliary", 3.9, 950),
]


@pytest.mark.parametrize(
    ("name", "pos", "ref", "alt", "edits", "rule", "subgroup", "percent", "count"),
    CASES,
    ids=[case[0] for case in CASES],
)
def test_each_terminal(name, pos, ref, alt, edits, rule, subgroup, percent, count):
    _expect(_call(pos, ref, alt, edits), rule, subgroup, percent, count)


def test_mid_exon_and_indel_stay_blank():
    missense = _call(300, "A", "G")
    assert pd.isna(missense["splice_rule"])
    assert pd.isna(missense["splice_sav_pct"])
    indel = pd.DataFrame(
        [{"chrom": "1", "pos": 300, "ref": "AT", "alt": "A", "gene": "GENE", "hgvsc": "", "evidence_score": 1}]
    )
    out = annotate_splice_heuristics(indel, [PLUS], MemoryFasta({"1": "A" * 800})).iloc[0]
    assert pd.isna(out["splice_rule"])
    assert out["evidence_score"] == 1


def test_minus_strand_donor_uses_coding_strand_ref():
    # Coding donor GT at +1/+2 is genomic CA (FASTA/VCF are plus-strand).
    # Coding G>A at +1 is genomic C>T.
    row = _call(600, "C", "T", {599: "A"}, exons=[MINUS])
    _expect(row, "DD2", "standard", 99.9, 2218)


def test_minus_strand_plus5_ref_g_is_genomic_c():
    # Coding +5 G (DD4 when E-1 is not G) is genomic C at exon.start-4 = 596.
    # Coding G>A is genomic C>T. Default E-1 genomic A -> coding T (not G).
    row = _call(596, "C", "T", exons=[MINUS])
    _expect(row, "DD4", "standard", 99.7, 583)


def test_minus_strand_genomic_gt_is_not_treated_as_coding_donor():
    # Genomic GT at +1/+2 is coding CA — not a GT donor, so no DD2 call.
    row = _call(600, "G", "A", {599: "T"}, exons=[MINUS])
    assert pd.isna(row["splice_rule"])


def test_gtf_exon_and_wrapped_fasta(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        'chr1\thavana\texon\t201\t400\t.\t+\t.\tgene_name "GBA"; transcript_id "ENST0001.5";\n',
        encoding="utf-8",
    )
    exons = load_gtf_exons(gtf, {"GBA1"}, {"GBA": "GBA1"})
    assert len(exons) == 1
    assert exons[0].start == 200 and exons[0].end == 400 and exons[0].gene == "GBA1"
    fasta = tmp_path / "ref.fa"
    # 406 bases, 8 per line, so +1/+2 (401-402) cross a line break. Header ">1\n" is 3 bytes.
    sequence = "A" * 400 + "GT" + "A" * 4
    lines = [sequence[index : index + 8] for index in range(0, len(sequence), 8)]
    fasta.write_text(">1\n" + "\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    index = tmp_path / "ref.fa.fai"
    index.write_text(f"1\t{len(sequence)}\t3\t8\t9\n", encoding="utf-8")
    reader = IndexedFasta(fasta, index)
    try:
        assert reader.fetch("chr1", 401, 402) == "GT"
        row = annotate_splice_heuristics(
            pd.DataFrame(
                [
                    {
                        "chrom": "1",
                        "pos": 401,
                        "ref": "G",
                        "alt": "A",
                        "gene": "GBA1",
                        "hgvsc": "ENST0001.5:c.100+1G>A",
                        "evidence_score": 4,
                    }
                ]
            ),
            exons,
            reader,
        ).iloc[0]
    finally:
        reader.close()
    _expect(row, "DD2", "standard", 99.9, 2218)


def test_missing_reference_leaves_columns_blank():
    frame = pd.DataFrame([{"chrom": "1", "pos": 401, "ref": "G", "alt": "A", "gene": "GENE"}])
    out = annotate_splice_heuristics(frame, None, None)
    assert pd.isna(out.loc[0, "splice_rule"])
    assert "splice_heuristic_n" in out.columns
