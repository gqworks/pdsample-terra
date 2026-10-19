"""AVI join and the local AVI table lookup script."""

from pathlib import Path

import pandas as pd
import pytest

from pdvar.assess import assess_sample
from pdvar.avi import merge_avi
from pdvar.config import load_config
import annotate_avi


def _vcf(path: Path, rows) -> Path:
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines += [f"{c}\t{p}\t.\t{r}\t{a}\t.\t{f}\t." for c, p, r, a, f in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


TABLE_ROWS = [
    ("chr1", 10, "A", "C", "0.10", "5.0"),
    ("chr1", 10, "A", "G", "1.25", "20.0"),
    ("chr1", 10, "A", "T", "0.20", "6.0"),
    ("chr1", 12, "C", "A", "0.30", "7.0"),
    ("chr1", 12, "C", "G", "2.50", "30.0"),
    ("chr1", 12, "C", "T", ".", "."),
    ("chr1", 5000, "G", "A", "0.40", "12.5"),
    ("chrM", 14, "G", "A", "0.05", "1.0"),
]


def _table(path: Path, rows=TABLE_ROWS, header="#CHROM\tPOS\tREF\tALT\traw_score\tPHRED", bare=False) -> Path:
    """bgzip + Tabix a small AVI table the way the download is laid out."""
    import pysam

    plain = path.with_suffix("")
    lines = [header] if header else []
    for chrom, pos, ref, alt, raw, phred in rows:
        name = chrom[3:] if bare else chrom
        name = "MT" if bare and name == "M" else name
        lines.append(f"{name}\t{pos}\t{ref}\t{alt}\t{raw}\t{phred}")
    plain.write_text("\n".join(lines) + "\n", encoding="utf-8")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
    return path


def test_phred_to_quantile_inverts_published_tail_ranks():
    assert annotate_avi.phred_to_quantile(20) == 0.99
    assert annotate_avi.phred_to_quantile(30) == 0.999
    assert annotate_avi.phred_to_quantile(0) == 0.0


def test_atlas_names_are_chr_prefixed():
    assert annotate_avi.atlas_chrom("1") == "chr1"
    assert annotate_avi.atlas_chrom("chrM") == "chrM"
    assert annotate_avi.atlas_chrom("MT") == "chrM"
    assert "q=chr1:10:A>G" in annotate_avi.atlas_url("1", 10, "A", "G")


def test_iter_snvs_keeps_pass_snvs_only(tmp_path: Path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t10\t.\tA\tG\t.\tPASS\t.",
                "chr1\t11\t.\tA\tGC\t.\tPASS\t.",
                "chr1\t12\t.\tC\tT,G\t.\tPASS\t.",
                "chr2\t13\t.\tA\tC\t.\tLowQual\t.",
                "MT\t14\t.\tG\tA\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    assert list(annotate_avi.iter_snvs(vcf)) == [
        ("chr1", 10, "A", "G"),
        ("chr1", 12, "C", "T"),
        ("chr1", 12, "C", "G"),
        ("MT", 14, "G", "A"),
    ]


def test_no_table_writes_header_only(tmp_path: Path):
    vcf = _vcf(tmp_path / "sample.vcf", [("chr1", 10, "A", "G", "PASS")])
    out = tmp_path / "avi.tsv"
    assert annotate_avi.main(["--vcf", str(vcf), "--out", str(out)]) == 0
    text = out.read_text(encoding="utf-8")
    assert text.startswith("chrom\tpos\tref\talt\tavi_raw\tavi_quantile\tavi_phred\tavi_url\n")
    assert len(text.splitlines()) == 1


def test_lookup_matches_exact_alleles_in_vcf_order(tmp_path: Path):
    table = _table(tmp_path / "avi.tsv.gz")
    vcf = _vcf(
        tmp_path / "sample.vcf",
        [
            ("1", 12, "C", "G,T", "PASS"),  # T has no score in the table
            ("1", 10, "A", "G", "PASS"),
            ("1", 11, "A", "GC", "PASS"),  # indel, never looked up
            ("1", 5000, "G", "A", "LowQual"),  # not PASS
            ("1", 7000, "T", "C", "PASS"),  # absent from the table
            ("2", 10, "A", "G", "PASS"),  # contig absent from the table
            ("MT", 14, "G", "A", "PASS"),
        ],
    )
    rows = annotate_avi.lookup_snvs(vcf, table)
    assert [(r["chrom"], r["pos"], r["ref"], r["alt"]) for r in rows] == [
        ("chr1", 12, "C", "G"),
        ("chr1", 10, "A", "G"),
        ("chrM", 14, "G", "A"),
    ]
    assert rows[0]["avi_raw"] == 2.5 and rows[0]["avi_phred"] == 30.0
    assert rows[0]["avi_quantile"] == 0.999
    assert rows[1]["avi_raw"] == 1.25 and rows[1]["avi_quantile"] == 0.99
    assert "q=chr1:10:A>G" in rows[1]["avi_url"]


def test_merge_gap_does_not_change_results(tmp_path: Path):
    table = _table(tmp_path / "avi.tsv.gz")
    vcf = _vcf(
        tmp_path / "sample.vcf",
        [("chr1", 10, "A", "G", "PASS"), ("chr1", 12, "C", "G", "PASS"), ("chr1", 5000, "G", "A", "PASS")],
    )
    tight = annotate_avi.lookup_snvs(vcf, table, merge_gap=0)
    wide = annotate_avi.lookup_snvs(vcf, table, merge_gap=100_000)
    assert tight == wide
    assert len(tight) == 3


def test_bare_contig_table_and_separate_index(tmp_path: Path):
    table_dir, index_dir = tmp_path / "table", tmp_path / "index"
    table_dir.mkdir()
    index_dir.mkdir()
    table = _table(table_dir / "avi.tsv.gz", bare=True)
    index = index_dir / "avi.tsv.gz.tbi"
    (table_dir / "avi.tsv.gz.tbi").rename(index)
    vcf = _vcf(tmp_path / "sample.vcf", [("chr1", 10, "A", "G", "PASS"), ("chrM", 14, "G", "A", "PASS")])
    out = tmp_path / "avi.tsv"
    args = ["--vcf", str(vcf), "--out", str(out), "--table", str(table), "--index", str(index)]
    assert annotate_avi.main(args) == 0
    written = pd.read_csv(out, sep="\t")
    assert list(written["chrom"]) == ["chr1", "chrM"]
    assert list(written["avi_phred"]) == [20.0, 1.0]


def test_header_aliases_and_missing_columns(tmp_path: Path):
    aliased = _table(tmp_path / "a.tsv.gz", header="#chrom\tpos\tref\talt\tavi_score\tphred_score")
    vcf = _vcf(tmp_path / "sample.vcf", [("chr1", 10, "A", "G", "PASS")])
    assert annotate_avi.lookup_snvs(vcf, aliased)[0]["avi_phred"] == 20.0
    broken = _table(tmp_path / "b.tsv.gz", header="#CHROM\tPOS\tREF\tALT\traw_score\tcadd")
    with pytest.raises(ValueError, match="phred"):
        annotate_avi.lookup_snvs(vcf, broken)


def test_ref_disagreement_fails_as_build_mismatch(tmp_path: Path):
    rows = [("chr1", pos, "A", "G", "1.0", "20.0") for pos in range(100, 160)]
    table = _table(tmp_path / "avi.tsv.gz", rows=rows)
    vcf = _vcf(tmp_path / "sample.vcf", [("chr1", pos, "C", "G", "PASS") for pos in range(100, 160)])
    with pytest.raises(ValueError, match="GRCh38"):
        annotate_avi.lookup_snvs(vcf, table)


def test_missing_table_file_is_an_error(tmp_path: Path):
    vcf = _vcf(tmp_path / "sample.vcf", [("chr1", 10, "A", "G", "PASS")])
    args = ["--vcf", str(vcf), "--out", str(tmp_path / "o.tsv"), "--table", str(tmp_path / "nope.tsv.gz")]
    assert annotate_avi.main(args) == 1


def test_merge_avi_joins_on_bare_chrom():
    variants = pd.DataFrame(
        [{"chrom": "1", "pos": 10, "ref": "A", "alt": "G", "variant_id": "1-10-A-G", "consequence": "missense_variant"}]
    )
    avi = pd.DataFrame(
        [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "avi_raw": 1.2, "avi_quantile": 0.999, "avi_phred": 30, "avi_url": "https://example.test"}]
    )
    merged = merge_avi(variants, avi)
    assert merged.iloc[0]["avi_phred"] == 30
    assert merged.iloc[0]["variant_id"] == "1-10-A-G"
    empty = merge_avi(variants, pd.DataFrame())
    assert pd.isna(empty.iloc[0]["avi_phred"])


def test_assess_uses_avi_as_evidence(tmp_path: Path):
    from test_assess import PANEL, _bed, _eh, _mt, _sv, _vep

    vep = tmp_path / "vep.vcf"
    _vep(vep)
    eh, mt, sv = tmp_path / "eh.vcf", tmp_path / "mt.vcf", tmp_path / "sv.vcf"
    _eh(eh)
    _mt(mt)
    _sv(sv)
    gene_bed, exon_bed = tmp_path / "genes.bed", tmp_path / "exons.bed"
    _bed(gene_bed, "PINK1")
    _bed(exon_bed, "PINK1")
    gba = tmp_path / "gauchian.json"
    gba.write_text('{"sample_id": "SYNTHETIC_SAMPLE", "GBA_copy_number": 0}\n', encoding="utf-8")
    introme = tmp_path / "introme.tsv"
    introme.write_text("CHROM\tPOS\tREF\tALT\tgene\tsamples_with_alt\tintrome_score\n", encoding="utf-8")
    avi = tmp_path / "avi.tsv"
    avi.write_text(
        "chrom\tpos\tref\talt\tavi_raw\tavi_quantile\tavi_phred\tavi_url\n"
        "chr1\t20650000\tG\tA\t2.5\t0.9995\t33\thttps://deepmind.google.com/science/alphagenome/atlas?m=variant&q=chr1:20650000:G>A\n",
        encoding="utf-8",
    )
    sample = {
        "sample_id": "SYNTHETIC_SAMPLE",
        "sex": "female",
        "vep_vcf": str(vep),
        "sv_vcf": str(sv),
        "expansionhunter_vcf": str(eh),
        "gauchian_json": str(gba),
        "mtdna_vcf": str(mt),
        "introme_tsv": str(introme),
        "avi_tsv": str(avi),
        "qc": {"pass": True},
    }
    written = assess_sample(
        sample,
        load_config(),
        panel_csv=PANEL,
        outdir=tmp_path / "out",
        gene_bed=gene_bed,
        exon_bed=exon_bed,
    )
    candidates = pd.read_csv(written["candidates_tsv"], sep="\t")
    hit = candidates[candidates["variant_id"] == "1-20650000-G-A"].iloc[0]
    assert "avi_high" in str(hit["evidence_reasons"])
    assert float(hit["avi_phred"]) == 33
    assert "alphagenome/atlas" in str(hit["avi_url"])
    report = written["report_md"].read_text(encoding="utf-8")
    assert "avi: called" in report
    assert "1-20650000-G-A" in report
