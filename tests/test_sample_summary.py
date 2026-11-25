"""APOE, m.2158T>C, Haplogrep, panel PRS, mean MLC, and MSS sample-summary fields."""

from pathlib import Path

import pytest

from pdvar.prs import load_prs_weights, score_prs
from pdvar.sample_summary import (
    apoe_genotype,
    build_sample_summary,
    mean_mlc,
    mt_2158_genotype,
    parse_haplogrep,
)
from pdvar.vcf_stream import chrom_bare

ROOT = Path(__file__).resolve().parents[1]


def _vcf(path: Path, chrom: str, pos: int, ref: str, alt: str, gt: str, extra_info: str = ".", fmt: str = "GT", sample: str = "S1") -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>',
                '##FORMAT=<ID=AF,Number=A,Type=Float,Description=AF>',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>',
                '##INFO=<ID=MI_MLC_score,Number=1,Type=Float,Description=MLC>',
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{extra_info}\t{fmt}\t{gt}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_apoe_e3e4(tmp_path):
    vcf = tmp_path / "apoe.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "chr19\t44908684\trs429358\tT\tC\t.\tPASS\t.\tGT\t0/1",
                "chr19\t44908822\trs7412\tC\tT\t.\tPASS\t.\tGT\t0/0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    from pdvar.prs import index_genotypes

    result = apoe_genotype(index_genotypes(vcf), "S1")
    assert result["apoe_genotype"] == "e3/e4"
    assert result["apoe_status"] == "called"


def test_apoe_unphased_e2e4_or_e1e3(tmp_path):
    vcf = tmp_path / "apoe.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "19\t44908684\t.\tT\tC\t.\tPASS\t.\tGT\t0/1",
                "19\t44908822\t.\tC\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    from pdvar.prs import index_genotypes

    result = apoe_genotype(index_genotypes(vcf), "S1")
    assert result["apoe_genotype"] == "e2/e4_or_e1/e3"


def test_mt_2158_heteroplasmy(tmp_path):
    vcf = tmp_path / "mt.vcf"
    _vcf(vcf, "chrM", 2158, "T", "C", "0/1:0.4:200", extra_info=".", fmt="GT:AF:DP")
    result = mt_2158_genotype(vcf, "S1", True)
    assert result["mt_2158_genotype"] == "T/C"
    assert result["mt_2158_heteroplasmy"] == 0.4


def test_mt_2158_absent_is_reference(tmp_path):
    vcf = tmp_path / "mt.vcf"
    _vcf(vcf, "chrM", 11778, "G", "A", "0/1:0.9:100", extra_info="MI_MLC_score=0.7", fmt="GT:AF:DP")
    result = mt_2158_genotype(vcf, "S1", True)
    assert result["mt_2158_genotype"] == "T"
    assert result["mt_2158_heteroplasmy"] == 0.0


def test_parse_haplogrep_found_polys(tmp_path):
    path = tmp_path / "haplogrep.tsv"
    path.write_text(
        "SampleID\tRange\tHaplogroup\tQuality\tFoundPolys\tNotFoundPolys\tRemainingPolys\n"
        "S1\t1-16569\tJ1b\t0.95\t73G 2158C 11719A\t\t1234C\n",
        encoding="utf-8",
    )
    parsed = parse_haplogrep(path)
    assert parsed["haplogroup"] == "J1b"
    assert parsed["haplogroup_quality"] == 0.95
    assert parsed["haplotag_positions"] == [73, 2158, 11719]
    assert parsed["haplotag_source"] == "haplogrep_found_polys"


def test_mean_mlc_excludes_haplotag(tmp_path):
    vcf = tmp_path / "mt.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=MI_MLC_score,Number=1,Type=Float,Description=MLC>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "chrM\t73\t.\tA\tG\t.\tPASS\tMI_MLC_score=0.1\tGT\t1",
                "chrM\t1234\t.\tC\tT\t.\tPASS\tMI_MLC_score=0.9\tGT\t1",
                "chrM\t2158\t.\tT\tC\t.\tPASS\tMI_MLC_score=0.5\tGT\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    result = mean_mlc(vcf, [73, 2158], haplotag_applied=True)
    assert result["mlc_n_excluded_haplotag"] == 2
    assert result["mean_mlc"] == 0.9
    assert result["mss"] == 0.9
    assert result["mlc_haplotag_exclusion"] == "haplogrep_found_polys"
    other = mean_mlc(vcf, [1234], haplotag_applied=True)
    assert other["mlc_n_excluded_haplotag"] == 1
    assert other["mean_mlc"] == 0.3
    assert other["mss"] == 0.6
    unused = mean_mlc(vcf, [73], haplotag_applied=False)
    assert unused["mlc_haplotag_exclusion"] == "not_applied"
    assert unused["mlc_n_excluded_haplotag"] == 0
    empty = mean_mlc(vcf, [], haplotag_applied=True)
    assert empty["mlc_haplotag_exclusion"] == "empty_found_polys"
    assert empty["mlc_n_excluded_haplotag"] == 0


def test_apoe_not_on_mendelian_panel():
    panel = (ROOT / "resources" / "panel" / "monopdaus-gene-inheritance.csv").read_text(encoding="utf-8")
    bed_rows = [
        line
        for line in (ROOT / "workflow" / "assets" / "summary_loci.grch38.bed").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    bed_text = "\n".join(bed_rows)
    assert "APOE" not in panel.upper()
    assert "FGF14" not in bed_text
    assert "NOTCH2NLC" not in bed_text
    assert "rs429358" in bed_text
    assert "rs7412" in bed_text


def test_load_prs_weights_rejects_uncited_rows(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text(
        "chrom\tpos\teffect_allele\tother_allele\tweight\trsid\tcitation\n"
        "12\t40340400\tA\tG\t2.0\trs34637584\t\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="citation"):
        load_prs_weights(path)


def test_panel_prs_weights_are_cited():
    weights = ROOT / "resources" / "prs" / "pgs000902_nalls2019.tsv"
    full = load_prs_weights(weights)
    panel = load_prs_weights(weights, panel_only=True)
    assert full is not None and panel is not None
    assert len(full) == 90
    assert len(panel) == 11
    assert full["citation"].astype(str).str.contains("Nalls").all()
    assert "GBA1" in set(panel["panel_gene"].astype(str))
    assert "SNCA" in set(panel["panel_gene"].astype(str))
    assert "LRRK2" in set(panel["panel_gene"].astype(str))
    bed_rsids = {
        line.split("\t")[3]
        for line in (ROOT / "workflow" / "assets" / "summary_loci.grch38.bed").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    }
    assert set(full["rsid"].astype(str)) <= bed_rsids
    assert {"rs429358", "rs7412"} <= bed_rsids
    assert len(bed_rsids) == 92


def test_score_prs_dosage(tmp_path):
    weights = tmp_path / "w.tsv"
    weights.write_text(
        "chrom\tpos\teffect_allele\tother_allele\tweight\trsid\tcitation\tscore_in_panel_prs\n"
        "12\t40340400\tA\tG\t2.0\trs34637584\tNalls 2019\ttrue\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "s.vcf"
    _vcf(vcf, "chr12", 40340400, "G", "A", "0/1")
    from pdvar.prs import index_genotypes

    scored = score_prs(load_prs_weights(weights), index_genotypes(vcf), "S1")
    assert scored["panel_prs"] == 2.0
    assert scored["prs_sites_used"] == 1
    assert scored["prs_sites_missing"] == 0


def test_build_sample_summary_end_to_end(tmp_path):
    summary_vcf = tmp_path / "summary.vcf"
    summary_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYN",
                "chr19\t44908684\trs429358\tT\tC\t.\tPASS\t.\tGT\t0/0",
                "chr19\t44908822\trs7412\tC\tT\t.\tPASS\t.\tGT\t0/0",
                "chr12\t40340400\trs34637584\tG\tA\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mt = tmp_path / "mt.vcf"
    mt.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=MI_MLC_score,Number=1,Type=Float,Description=MLC>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=AF,Number=A,Type=Float,Description=AF>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYN",
                "chrM\t73\t.\tA\tG\t.\tPASS\tMI_MLC_score=0.2\tGT:AF\t1:1",
                "chrM\t2158\t.\tT\tC\t.\tPASS\tMI_MLC_score=0.4\tGT:AF\t0/1:0.3",
                "chrM\t1234\t.\tC\tT\t.\tPASS\tMI_MLC_score=0.8\tGT:AF\t1:1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    hap = tmp_path / "hap.tsv"
    hap.write_text(
        "SampleID\tHaplogroup\tQuality\tFoundPolys\nS1\tH1\t0.9\t73G\n",
        encoding="utf-8",
    )
    sample = {
        "sample_id": "SYN",
        "summary_vcf": str(summary_vcf),
        "mtdna_vcf": str(mt),
        "mtdna_annotation": "called",
        "haplogrep_tsv": str(hap),
        "prs_weights": str(ROOT / "resources" / "prs" / "pgs000902_nalls2019.tsv"),
    }
    summary = build_sample_summary(sample)
    assert summary["apoe_genotype"] == "e3/e3"
    assert summary["mt_2158_genotype"] == "T/C"
    assert summary["haplogroup"] == "H1"
    assert summary["panel_prs"] == 2.4289
    assert summary["prs_sites_total"] == 90
    assert summary["prs_sites_used"] == 1
    assert summary["mean_mlc"] == 0.6
    assert summary["mss"] == 1.2
    assert summary["mlc_n_excluded_haplotag"] == 1
    assert summary["mlc_haplotag_exclusion"] == "haplogrep_found_polys"
    assert chrom_bare("chrM") == "M"
