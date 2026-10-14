"""AVI join and the Atlas lookup script, without calling the API."""

import os
from pathlib import Path

import pandas as pd

from pdvar.assess import assess_sample
from pdvar.avi import merge_avi
from pdvar.config import load_config
import annotate_avi


def test_phred_matches_published_tail_ranks():
    assert annotate_avi.avi_phred(0.99) == 20
    assert annotate_avi.avi_phred(0.999) == 30
    assert annotate_avi.avi_phred(1.0) == 100


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


def test_empty_key_writes_header_only(tmp_path: Path, monkeypatch):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tG\t.\tPASS\t.\n", encoding="utf-8")
    out = tmp_path / "avi.tsv"
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    assert annotate_avi.main(["--vcf", str(vcf), "--out", str(out)]) == 0
    text = out.read_text(encoding="utf-8")
    assert text.startswith("chrom\tpos\tref\talt\tavi_raw\tavi_quantile\tavi_phred\tavi_url\n")
    assert len(text.splitlines()) == 1
    assert os.environ.get("ALPHAGENOME_API_KEY", "") == ""


class _Cell:
    def __init__(self, value):
        self.value = value

    def __getitem__(self, _index):
        return self.value


class _Obs:
    def __init__(self, variant: str):
        self._variant = variant

    def __getitem__(self, key):
        assert key == "variant"
        return type("Col", (), {"iloc": [self._variant]})()


class _Scores:
    def __init__(self, raw, cdf, variant):
        self.n_obs = 1
        self.X = _Cell(raw)
        self.layers = {"quantiles": _Cell(cdf)}
        self.obs = _Obs(variant)

    def get(self, key):
        return self if key == "AVI_SCORE" else None


def test_row_rejects_a_different_returned_allele():
    scores = {"AVI_SCORE": _Scores(1.0, 0.99, "chr1:10:A>T")}
    assert annotate_avi._row_from_scores("chr1", 10, "A", "G", scores) is None


def test_row_from_matching_score():
    scores = {"AVI_SCORE": _Scores(1.25, 0.99, "chr1:10:A>G")}
    row = annotate_avi._row_from_scores("1", 10, "A", "G", scores)
    assert row is not None
    assert row["avi_phred"] == 20
    assert row["avi_raw"] == 1.25
    assert row["chrom"] == "chr1"


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
