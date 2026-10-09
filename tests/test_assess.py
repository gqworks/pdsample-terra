"""End-to-end assessment on one synthetic sample."""

import json
from pathlib import Path

import pandas as pd
import pytest

from pdvar.assess import assess_sample
from pdvar.config import load_config

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "resources" / "panel" / "monopdaus-gene-inheritance.csv"

CSQ = (
    "Allele|Consequence|SYMBOL|HGVSc|HGVSp|CANONICAL|CADD_PHRED|REVEL|"
    "SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL"
)


def _vep(path: Path) -> None:
    header = [
        "##fileformat=VCFv4.2",
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {CSQ}">',
        '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD AF">',
        '##INFO=<ID=AF_popmax,Number=A,Type=Float,Description="popmax AF">',
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=GQ>",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=AD>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC",
    ]
    rows = [
        "1\t20650000\t.\tG\tA\t.\tPASS\t"
        "CSQ=A|missense_variant|PINK1|c.1G>A|p.Met1Ile|YES|30|0.9||||;gnomAD_AF=0;AF_popmax=0\t"
        "GT:GQ:DP:AD\t0/1:99:40:20,20",
        "1\t20651000\t.\tC\tT\t.\tPASS\t"
        "CSQ=T|stop_gained|PINK1|c.2C>T|p.Gln1Ter|YES|40|0||||;gnomAD_AF=0;AF_popmax=0\t"
        "GT:GQ:DP:AD\t0/1:99:40:20,20",
        "1\t20652000\t.\tA\tG\t.\tPASS\t"
        "CSQ=G|missense_variant|PINK1|c.3A>G|p.X|YES|25|0.8||||;gnomAD_AF=0.5;AF_popmax=0.5\t"
        "GT:GQ:DP:AD\t0/1:99:40:20,20",
        "2\t179390000\t.\tC\tT\t.\tPASS\t"
        "CSQ=T|stop_gained|TTN|c.10C>T|p.Gln4Ter|YES|40|||||;gnomAD_AF=0\t"
        "GT:GQ:DP:AD\t0/1:99:40:20,20",
    ]
    path.write_text("\n".join(header + rows) + "\n", encoding="utf-8")


def _eh(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=REPID,Number=1,Type=String,Description=Locus>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=REPCN,Number=1,Type=String,Description=Repeats>",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_SAMPLE",
                "4\t39348425\t.\tA\t<STR>\t.\tPASS\tREPID=RFC1\tGT:REPCN:DP\t0/1:11/12:30",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _mt(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=GENE,Number=1,Type=String,Description=Gene>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=AF,Number=A,Type=Float,Description=AF>",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_SAMPLE",
                "M\t11778\t.\tG\tA\t.\tPASS\tGENE=MT-ND4\tGT:AF:DP\t0/1:0.4:100",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _sv(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_SAMPLE",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _bed(path: Path, gene: str) -> None:
    path.write_text(f"chr1\t1000\t2000\t.\t{gene}\n", encoding="utf-8")


def test_assess_ranks_one_synthetic_sample(tmp_path):
    vep = tmp_path / "vep.vcf"
    eh = tmp_path / "eh.vcf"
    mt = tmp_path / "mt.vcf"
    sv = tmp_path / "sv.vcf"
    gene_bed = tmp_path / "genes.bed"
    exon_bed = tmp_path / "exons.bed"
    gba = tmp_path / "gauchian.json"
    introme = tmp_path / "introme.tsv"
    _vep(vep)
    _eh(eh)
    _mt(mt)
    _sv(sv)
    _bed(gene_bed, "PINK1")
    _bed(exon_bed, "PINK1")
    gba.write_text(
        json.dumps({"sample_id": "SYNTHETIC_SAMPLE", "GBA_copy_number": 0}) + "\n",
        encoding="utf-8",
    )
    introme.write_text("CHROM\tPOS\tREF\tALT\tgene\tsamples_with_alt\tintrome_score\n", encoding="utf-8")
    sample = {
        "sample_id": "SYNTHETIC_SAMPLE",
        "sex": "female",
        "vep_vcf": str(vep),
        "sv_vcf": str(sv),
        "expansionhunter_vcf": str(eh),
        "gauchian_json": str(gba),
        "mtdna_vcf": str(mt),
        "introme_tsv": str(introme),
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
    by_gene = candidates.groupby("gene")["inheritance_status"].apply(lambda values: set(values))
    assert by_gene["PINK1"] == {"AR_comphet_candidate"}
    assert "1-20652000-A-G" not in set(candidates["variant_id"].astype(str))
    assert by_gene["GBA1"] == {"risk_biallelic"}
    ttn = candidates[candidates["gene"] == "TTN"].iloc[0]
    assert ttn["inheritance_status"] == "incidental"
    assert ttn["incidental"] == "TTN_truncation"
    assert ttn["clinvar_class"] != "B/LB"
    assert by_gene["MT-ND4"] == {"mtdna_variant"}
    assert "RFC1" not in set(candidates["gene"].astype(str))

    liability = pd.read_csv(written["liability_tsv"], sep="\t")
    row = liability.iloc[0]
    assert row["classification"] == "solved"
    assert row["solved_genes"] == "PINK1"
    assert "GBA1" in str(row["partial_genes"])
    assert pd.isna(row["explained_liability"])

    qc = json.loads(written["qc_json"].read_text(encoding="utf-8"))
    assert qc["sv_status"] == "called"
    assert qc["modules"]["gatk_sv"] == "called"
    assert qc["sv"][0]["sv_status"] == "no SV"
    callsets = {item["module"]: item["callset_status"] for item in qc["callsets"]}
    assert callsets["expansionhunter"] == "STR carrier"
    assert "SYNTHETIC_SAMPLE" in written["report_md"].read_text(encoding="utf-8")
    assert written["report_html"].read_text(encoding="utf-8").startswith("<!DOCTYPE html>")
    assert written["followup_tsv"].exists()


def test_assess_requires_sv_vcf(tmp_path):
    vep = tmp_path / "vep.vcf"
    _vep(vep)
    sample = {"sample_id": "SYNTHETIC_SAMPLE", "sex": "female", "vep_vcf": str(vep)}
    with pytest.raises(FileNotFoundError, match="sv_vcf"):
        assess_sample(sample, load_config(), panel_csv=PANEL, outdir=tmp_path / "out")
