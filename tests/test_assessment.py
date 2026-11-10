"""SV dosage, liability, and biomarkers."""

from pathlib import Path

import pandas as pd
import pytest

from pdvar import biomarkers, liability, sv
from pdvar.config import load_config
from pdvar.splicing import merge_splicevault, summarize_splice_scores


def _beds(directory: Path) -> tuple[Path, Path]:
    genes = directory / "genes.bed"
    exons = directory / "exons.bed"
    genes.write_text("chr4\t800\t2500\tENSG\tSNCA\nchr6\t900\t4000\tENSG\tPRKN\n", encoding="utf-8")
    exons.write_text("chr4\t1000\t2000\tENSG\tSNCA\nchr6\t1000\t1200\tENSG\tPRKN\n", encoding="utf-8")
    return genes, exons


def test_depth_only_artifact_helpers():
    info = {"ALGORITHMS": "depth", "AF": "0.49", "SVTYPE": "DEL"}
    fmt = {"GT": "1/1", "EV": "RD", "RD_GQ": "1", "OGQ": "1", "GQ": "43"}
    assert sv.is_depth_only_algorithms(info)
    assert not sv.has_pesr_format_evidence(fmt)
    assert sv.is_depth_only_artifact(info, fmt, 0.49)
    assert not sv.is_depth_only_artifact(info, fmt, 0.01)  # rare depth CNV kept
    assert not sv.is_depth_only_artifact(info, {**fmt, "EV": "RD,PE"}, 0.49)
    assert not sv.is_depth_only_artifact(info, {**fmt, "RD_GQ": "20", "OGQ": "20"}, 0.49)
    assert not sv.is_depth_only_algorithms({"ALGORITHMS": "manta,depth"})
    assert sv.is_depth_only_artifact(info, {"EV": "RD", "OGQ": "2"}, 0.2)  # OGQ-only weak


def test_sv_dosage_and_not_assessed(tmp_path):
    genes, exons = _beds(tmp_path)
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Type>",
                "##INFO=<ID=PREDICTED_COPY_GAIN,Number=.,Type=String,Description=Gain>",
                "##INFO=<ID=PREDICTED_LOF,Number=.,Type=String,Description=LoF>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=CN>",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=ECN>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_C",
                "chr4\t900\t.\tN\t<DUP>\t.\tPASS\tEND=2200;SVTYPE=DUP;PREDICTED_COPY_GAIN=SNCA\tGT:CN:ECN\t0/0:2:2\t0/1:4:2",
                "chr6\t1000\tdel\tN\t<DEL>\t.\tPASS\tEND=1300;SVTYPE=DEL;PREDICTED_LOF=PRKN\tGT:CN:ECN\t./.:1:2\t0/0:2:2",
                "chr6\t2000\tmulti\tN\t<DEL>\t.\tMULTIALLELIC\tEND=2300;SVTYPE=DEL;PREDICTED_LOF=PRKN\tGT:CN:ECN\t0/1:1:2\t0/0:2:2",
                "chr6\t1000\tdrop\tN\t<DEL>\t.\tHIGH_NCR\tEND=1300;SVTYPE=DEL;PREDICTED_LOF=PRKN\tGT:CN:ECN\t0/1:1:2\t0/0:2:2",
                "chr6\t1000\tdepth\tN\t<DEL>\t.\tDEPTH_LT_5KB;SAMPLE_WGD_OUTLIER\tEND=1300;SVTYPE=DEL;PREDICTED_LOF=PRKN\tGT:CN:ECN\t0/1:1:2\t0/0:2:2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    panel = pd.DataFrame({"gene": ["SNCA", "PRKN"], "gene_id": ["1", "2"]})
    table, assessment = sv.ingest_sv(
        vcf, genes, exons, panel, snv_samples=["SAMPLE_A", "SAMPLE_C", "SAMPLE_D"], auto_index=False
    )
    snca = table[(table["gene"] == "SNCA") & (table["sample_id"] == "SAMPLE_C")].iloc[0]
    assert snca["dosage_event"] == "triplication_or_higher"
    assert snca["consequence"] == "copy_number_gain"
    assert snca["zygosity"] == "het"
    prkn = table[(table["gene"] == "PRKN") & (table["sample_id"] == "SAMPLE_A") & (table["sv_id"] == "del")].iloc[0]
    assert prkn["zygosity"] == "het"
    assert prkn["cn_inferred_gt"] == True  # noqa: E712
    assert prkn["consequence"] == "transcript_ablation"
    assert prkn["svlen"] == 300
    assert snca["svlen"] == 1300
    assert sv.sv_length({"SVLEN": "-250"}, 1000, 1250) == -250
    assert sv.sv_length({}, 1000, 1000) is None
    assert prkn["filter"] == "PASS"
    multi = table[(table["gene"] == "PRKN") & (table["sv_id"] == "multi")].iloc[0]
    assert multi["filter"] == "MULTIALLELIC"
    assert "del" in set(table["sv_id"].astype(str))
    assert "multi" in set(table["sv_id"].astype(str))
    assert "drop" not in set(table["sv_id"].astype(str))
    assert "depth" not in set(table["sv_id"].astype(str))
    assert set(table["filter"].astype(str)) <= {"PASS", "MULTIALLELIC"}
    status = assessment.set_index("sample_id")["sv_status"]
    assert status["SAMPLE_D"] == "SV not assessed"
    assert status["SAMPLE_C"] == "SV carrier"


def test_sv_drops_common_depth_only_artifacts(tmp_path):
    """SOX6-style PASS depth DELs (no PE/SR, weak RD_GQ, high cohort AF) are dropped."""
    genes = tmp_path / "genes.bed"
    exons = tmp_path / "exons.bed"
    genes.write_text(
        "chr11\t15956448\t16749591\tENSG\tSOX6\nchr6\t900\t4000\tENSG\tPRKN\n",
        encoding="utf-8",
    )
    exons.write_text(
        "chr11\t16231128\t16234771\tENSG\tSOX6\nchr6\t1000\t1200\tENSG\tPRKN\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Type>",
                "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=Len>",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=AF>",
                "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=Algs>",
                "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=Evidence>",
                "##INFO=<ID=PREDICTED_LOF,Number=.,Type=String,Description=LoF>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=ECN>",
                "##FORMAT=<ID=EV,Number=1,Type=String,Description=EV>",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=GQ>",
                "##FORMAT=<ID=OGQ,Number=1,Type=Integer,Description=OGQ>",
                "##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description=RD_CN>",
                "##FORMAT=<ID=RD_GQ,Number=1,Type=Integer,Description=RD_GQ>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A",
                # SOX6-like artifact: depth-only, EV=RD, RD_GQ=1, AF=0.49
                "chr11\t16228501\tsox6_art\tN\t<DEL>\t55\tPASS\tEND=16239001;SVTYPE=DEL;SVLEN=10501;AF=0.488532;ALGORITHMS=depth;EVIDENCE=BAF,RD;PREDICTED_LOF=SOX6\tGT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ\t1/1:2:RD:43:1:0:1",
                # Rare depth-only CNV kept (AF below threshold)
                "chr6\t1000\trare_depth\tN\t<DEL>\t.\tPASS\tEND=1300;SVTYPE=DEL;SVLEN=300;AF=0.002;ALGORITHMS=depth;EVIDENCE=RD;PREDICTED_LOF=PRKN\tGT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ\t0/1:2:RD:60:1:1:1",
                # Common but PE-supported — kept
                "chr6\t2000\tpesr\tN\t<DEL>\t.\tPASS\tEND=2300;SVTYPE=DEL;SVLEN=300;AF=0.4;ALGORITHMS=manta;EVIDENCE=PE,SR;PREDICTED_LOF=PRKN\tGT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ\t0/1:2:PE,SR:90:90:1:2",
                # Depth-only common but strong RD_GQ/OGQ — kept
                "chr6\t3000\tstrong_depth\tN\t<DEL>\t.\tPASS\tEND=3300;SVTYPE=DEL;SVLEN=300;AF=0.4;ALGORITHMS=depth;EVIDENCE=RD;PREDICTED_LOF=PRKN\tGT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ\t0/1:2:RD:80:20:1:20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    panel = pd.DataFrame({"gene": ["SOX6", "PRKN"], "gene_id": ["1", "2"]})
    table, _ = sv.ingest_sv(vcf, genes, exons, panel, snv_samples=["SAMPLE_A"], auto_index=False)
    ids = set(table["sv_id"].astype(str))
    assert "sox6_art" not in ids
    assert "rare_depth" in ids
    assert "pesr" in ids
    assert "strong_depth" in ids
    rare = table[table["sv_id"] == "rare_depth"].iloc[0]
    assert rare["sv_algorithms"] == "depth"
    assert rare["sv_evidence"] == "RD"
    assert float(rare["rd_gq"]) == 1.0


def test_liability_stays_null_without_citation():
    variants = pd.DataFrame(
        [
            {"sample_id": "SAMPLE_A", "gene": "PRKN", "inheritance_status": "AR_hom", "variant_id": "v", "variant_class": "SNV"},
            {"sample_id": "SAMPLE_B", "gene": "LRRK2", "inheritance_status": "AD_candidate", "variant_id": "v2", "variant_class": "SNV"},
        ]
    )
    plain = liability.classify_cohort(variants, ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"], cfg={"reduced_penetrance_genes": ["LRRK2"]})
    by_sample = plain.set_index("sample_id")
    assert by_sample.loc["SAMPLE_A", "classification"] == "solved"
    assert pd.isna(by_sample.loc["SAMPLE_A", "explained_liability"])
    assert by_sample.loc["SAMPLE_B", "classification"] == "partially_explained"
    assert by_sample.loc["SAMPLE_C", "classification"] == "unexplained"
    cited = liability.classify_cohort(
        variants.iloc[:1],
        ["SAMPLE_A"],
        cfg={"penetrance": {"PRKN": {"biallelic": 0.9, "citation": "Curator note, test only."}}},
    )
    assert cited.iloc[0]["explained_liability"] == pytest.approx(0.9)


def test_prs_weights_require_citation(tmp_path):
    path = tmp_path / "weights.tsv"
    path.write_text("chrom\tpos\teffect_allele\tweight\tcitation\n1\t10\tA\t0.2\t\n", encoding="utf-8")
    with pytest.raises(ValueError, match="citation"):
        liability.load_prs_weights(path)
    empty = tmp_path / "empty.tsv"
    empty.write_text("chrom\tpos\teffect_allele\tweight\tcitation\n", encoding="utf-8")
    assert liability.load_prs_weights(empty) is None
    shipped = load_config()
    assert shipped["liability"]["prs_beta"] is None


def test_biomarker_multiplier_needs_citation():
    df = pd.DataFrame(
        [{"sample_id": "SAMPLE_A", "gene": "PRKN", "age_onset": 30, "inheritance_status": "AR_hom", "phenotype_score": 4}]
    )
    cfg = {
        "gene_classes": {"recessive_eopd": ["PRKN"]},
        "fields": {"age_onset": ["age_onset"]},
        "rules": [
            {"id": "early", "field": "age_onset", "gene_class": "recessive_eopd", "op": "below", "value": 50, "multiplier": 2, "citation": None}
        ],
    }
    out = biomarkers.apply_biomarker_priors(df, cfg)
    assert out.iloc[0]["biomarker_prior_multiplier"] == 1
    assert out.iloc[0]["biomarker_prior_adjusted"] == False  # noqa: E712
    assert out.iloc[0]["phenotype_score"] == 4
    assert "missing_citation" in out.iloc[0]["biomarker_notes"]


def test_splicevault_keeps_highest_proportion_without_new_points():
    variants = pd.DataFrame(
        [{"variant_id": "1-10-A-G", "evidence_score": 4, "evidence_reasons": "introme_high", "consequence": "splice_region_variant"}]
    )
    vault = pd.DataFrame(
        [
            {"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "gene": "PINK1", "event": "skip", "outcome": "exon_skip", "proportion": 0.2, "exon_start": 1, "exon_end": 5},
            {"chrom": "1", "pos": 10, "ref": "A", "alt": "G", "gene": "PINK1", "event": "cryptic", "outcome": "cryptic_donor", "proportion": 0.8, "exon_start": 8, "exon_end": 12},
        ]
    )
    merged = merge_splicevault(variants, vault)
    assert merged.iloc[0]["splicevault_outcome"] == "cryptic_donor"
    assert merged.iloc[0]["evidence_score"] == 4
    summary = summarize_splice_scores(pd.DataFrame([{"introme_pangolin_gain": 0.1, "introme_pangolin_loss": 0.4}]))
    assert summary.iloc[0]["pangolin_max"] == 0.4
    assert "alphamissense" in summary.columns


def test_large_vcf_is_not_scanned(tmp_path, monkeypatch):
    from pdvar import vcf_stream

    vcf = tmp_path / "big.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\nchr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    monkeypatch.setattr(vcf_stream, "STREAM_MAX_BYTES", 1)
    with pytest.raises(RuntimeError, match="Refusing to scan"):
        list(vcf_stream.iter_vcf(vcf, allow_full_scan=False))
    records = list(vcf_stream.iter_vcf(vcf, allow_full_scan=True))
    assert records[0]["pos"] == 1
