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
    genes.write_text("chr4\t800\t25000\tENSG\tSNCA\nchr6\t900\t4000\tENSG\tPRKN\n", encoding="utf-8")
    exons.write_text("chr4\t1000\t2000\tENSG\tSNCA\nchr6\t1000\t1200\tENSG\tPRKN\n", encoding="utf-8")
    return genes, exons


def test_depth_only_artifact_helpers():
    info = {"ALGORITHMS": "depth", "AF": "0.49", "SVTYPE": "DEL"}
    fmt = {"GT": "1/1", "EV": "RD", "RD_GQ": "1", "OGQ": "1", "GQ": "43"}
    assert sv.is_depth_only_algorithms(info)
    assert not sv.has_pesr_format_evidence(fmt)
    assert sv.is_depth_only_artifact(info, fmt)  # cohort AF is not required
    assert not sv.is_depth_only_artifact(info, {**fmt, "EV": "RD,PE"})
    assert not sv.is_depth_only_artifact(info, {**fmt, "RD_GQ": "20", "OGQ": "20"})
    assert not sv.is_depth_only_algorithms({"ALGORITHMS": "manta,depth"})
    assert sv.is_depth_only_artifact(info, {"EV": "RD", "OGQ": "2"})  # OGQ-only weak


def test_pesr_caller_artifact_helpers():
    """Raw Manta-like DELs without CN are dropped; GATK-SV CN records are kept."""
    manta = {"GT": "0/1", "GQ": "144", "PR": "14,4", "SR": "14,9"}
    assert sv.lacks_depth_cn(manta)
    assert not sv.lacks_depth_cn({"GT": "0/1", "CN": "1", "ECN": "2"})
    assert not sv.lacks_depth_cn({"GT": "0/1", "RD_CN": "1", "ECN": "2"})
    assert sv.depth_copy_number({"GT": "0/1", "RD_CN": "1"}) == 1
    assert sv.depth_copy_number({"CN": "2", "RD_CN": "1"}) == 2
    assert sv.manta_alt_count(manta, "SR") == 9
    assert sv.manta_alt_count({"SR": "3,0"}, "SR") == 0
    assert sv.manta_alt_count({"PR": "14,5"}, "SR") is None
    assert sv.event_span_bp({"SVLEN": "-97643468"}, 20604914, 118248382) == 97643468
    # 97 Mb spanning DEL
    assert sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL"}, manta, svtype="DEL", span_bp=97_643_468, zygosity="het"
    )
    # IMPRECISE + low JUNCTION_QUAL
    assert sv.is_pesr_caller_artifact(
        {"IMPRECISE": True, "JUNCTION_QUAL": "1"},
        {"GT": "0/1", "PR": "14,5"},
        svtype="DEL",
        span_bp=22_887_154,
        zygosity="het",
    )
    # Large DEL without SR
    assert sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL"},
        {"GT": "0/1", "PR": "14,5"},
        svtype="DEL",
        span_bp=100_000,
        zygosity="het",
    )
    # 1.66 Mb TWNK-like: SR-supported but no CN — require_cn_min_span_bp drops it
    assert sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL"},
        {"GT": "0/1", "GQ": "310", "PR": "3,4", "SR": "30,27"},
        svtype="DEL",
        span_bp=1_657_725,
        zygosity="het",
    )
    # Just under CN floor with SR — kept (e.g. 99 kb)
    assert not sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL"},
        {"GT": "0/1", "SR": "5,8", "PR": "4,3"},
        svtype="DEL",
        span_bp=99_000,
        zygosity="het",
    )
    # Homozygous multi-kb without CN
    assert sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL"},
        {"GT": "1/1", "SR": "0,6", "PR": "0,3"},
        svtype="DEL",
        span_bp=100_000,
        zygosity="hom",
    )
    # Small precise DEL with SR — kept
    assert not sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL", "CIPOS": "0,16"},
        {"GT": "0/1", "SR": "0,6", "PR": "0,3"},
        svtype="DEL",
        span_bp=800,
        zygosity="het",
    )
    # GATK-SV with CN skips PE/SR-caller QC even for multi-Mb span
    assert not sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL"},
        {"GT": "0/1", "CN": "1", "ECN": "2"},
        svtype="DEL",
        span_bp=5_000_000,
        zygosity="het",
    )
    # GATK-SV joint VCFs write RD_CN, not CN. 54 kb PRKN DEL must skip Manta PESR QC.
    assert not sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL", "ALGORITHMS": "depth,manta"},
        {"GT": "0/1", "RD_CN": "1", "ECN": "2", "EV": "RD,SR", "SR_GT": "1", "SR_GQ": "67"},
        svtype="DEL",
        span_bp=54_475,
        zygosity="het",
    )
    # Same span without CN or RD_CN still drops (no Manta SR=ref,alt)
    assert sv.is_pesr_caller_artifact(
        {"SVTYPE": "DEL"},
        {"GT": "0/1", "EV": "RD,SR", "SR_GT": "1"},
        svtype="DEL",
        span_bp=54_475,
        zygosity="het",
    )


def test_sv_drops_raw_manta_spanning_dels(tmp_path):
    """PASS Manta multi-Mb DELs without CN must not become transcript_ablation solves."""
    genes = tmp_path / "genes.bed"
    exons = tmp_path / "exons.bed"
    genes.write_text(
        "chr10\t10000000\t120000000\tENSG\tTWNK\n"
        "chr6\t158127000\t158129000\tENSG\tSERAC1\n"
        "chr6\t900\t4000\tENSG\tPRKN\n",
        encoding="utf-8",
    )
    exons.write_text(
        "chr10\t102839000\t102840000\tENSG\tTWNK\n"
        "chr6\t158127400\t158127800\tENSG\tSERAC1\n"
        "chr6\t1000\t1200\tENSG\tPRKN\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "manta.sv.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Type>",
                "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=Len>",
                "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=Imprecise>",
                "##INFO=<ID=JUNCTION_QUAL,Number=1,Type=Integer,Description=JQ>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=GQ>",
                "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=PR>",
                "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=SR>",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=CN>",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=ECN>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A",
                # APGS-like 97 Mb spanning DEL — drop
                "chr10\t20604914\thuge\tG\t<DEL>\t144\tPASS\tEND=118248382;SVTYPE=DEL;SVLEN=-97643468\tGT:GQ:PR:SR\t0/1:144:14,4:14,9",
                # IMPRECISE large DEL, JUNCTION_QUAL=1 — drop
                "chr10\t50000000\timprecise\tN\t<DEL>\t999\tPASS\tEND=80000000;SVTYPE=DEL;SVLEN=-30000000;IMPRECISE;JUNCTION_QUAL=1\tGT:GQ:PR\t0/1:450:14,5",
                # Hom 15 Mb — drop
                "chr10\t40000000\thom_huge\tN\t<DEL>\t940\tPASS\tEND=55000000;SVTYPE=DEL;SVLEN=-15000000\tGT:GQ:PR:SR\t1/1:29:0,3:3,23",
                # APGS_004451-like 1.66 Mb TWNK DEL with SR but no CN — drop
                "chr10\t99838115\ttwnk_mid\tN\t<DEL>\t999\tPASS\tEND=101495840;SVTYPE=DEL;SVLEN=-1657725\tGT:GQ:PR:SR\t0/1:310:3,4:30,27",
                # Small SERAC1-like 800 bp with SR — keep
                "chr6\t158127234\tsmall\tT\t<DEL>\t316\tPASS\tEND=158128034;SVTYPE=DEL;SVLEN=-800\tGT:GQ:PR:SR\t0/1:24:0,3:0,6",
                # GATK-SV-like with CN, multi-Mb — keep (CN path skips PE/SR QC)
                "chr6\t1000\tgatk_cn\tN\t<DEL>\t.\tPASS\tEND=3001000;SVTYPE=DEL;SVLEN=-3000000\tGT:CN:ECN\t0/1:1:2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    panel = pd.DataFrame({"gene": ["TWNK", "SERAC1", "PRKN"], "gene_id": ["1", "2", "3"]})
    table, _ = sv.ingest_sv(vcf, genes, exons, panel, snv_samples=["SAMPLE_A"], auto_index=False)
    ids = set(table["sv_id"].astype(str))
    assert "huge" not in ids
    assert "imprecise" not in ids
    assert "hom_huge" not in ids
    assert "twnk_mid" not in ids
    assert "small" in ids
    assert "gatk_cn" in ids
    small = table[table["sv_id"] == "small"].iloc[0]
    assert small["gene"] == "SERAC1"
    assert small["consequence"] == "transcript_ablation"


def test_gatk_sv_rd_cn_prkn_del_is_kept():
    """Joint GATK-SV PRKN exon DEL with RD_CN (no CN) must not be PESR-dropped."""
    from pdvar.config import PROJECT_ROOT

    vcf = Path(__file__).resolve().parent / "fixtures" / "sv" / "apgs_000502_prkn_del.vcf"
    genes = PROJECT_ROOT / "resources" / "panel" / "monopdaus-genes.grch38.bed"
    exons = PROJECT_ROOT / "resources" / "panel" / "monopdaus-exons.grch38.bed"
    panel = pd.DataFrame({"gene": ["PRKN"], "gene_id": ["ENSG00000185345"]})
    table, _ = sv.ingest_sv(vcf, genes, exons, panel, snv_samples=["APGS_000502"], auto_index=False)
    row = table[table["sv_id"].astype(str).str.contains("DEL_chr6_1438")].iloc[0]
    assert row["gene"] == "PRKN"
    assert int(row["copy_number"]) == 1
    assert int(row["exons_overlapped"]) >= 1
    assert row["consequence"] == "transcript_ablation"
    assert bool(row["scored_dosage"]) is True
    assert bool(row["sv_artifact"]) is False


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
                "chr4\t900\t.\tN\t<DUP>\t.\tPASS\tEND=20000;SVTYPE=DUP;PREDICTED_COPY_GAIN=SNCA\tGT:CN:ECN\t0/0:2:2\t0/1:4:2",
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
    assert snca["svlen"] == 19100
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
    """PASS depth DELs with no PE/SR and weak RD_GQ/OGQ stay as unscored artifacts."""
    genes = tmp_path / "genes.bed"
    exons = tmp_path / "exons.bed"
    genes.write_text(
        "chr11\t15956448\t16749591\tENSG\tSOX6\nchr6\t900\t4000\tENSG\tPRKN\nchr3\t193593000\t193690000\tENSG\tOPA1\n",
        encoding="utf-8",
    )
    exons.write_text(
        "chr11\t16231128\t16234771\tENSG\tSOX6\nchr6\t1000\t1200\tENSG\tPRKN\nchr3\t193663600\t193673500\tENSG\tOPA1\n",
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
                # OPA1-like rare depth-only with weak RD_GQ — also dropped (no AF gate)
                "chr3\t193663575\topa1_art\tN\t<DEL>\t.\tPASS\tEND=193673575;SVTYPE=DEL;SVLEN=10001;AF=0.020642;ALGORITHMS=depth;EVIDENCE=RD;PREDICTED_LOF=OPA1\tGT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ\t0/1:2:RD:58:2:1:1",
                # PE-supported — kept
                "chr6\t2000\tpesr\tN\t<DEL>\t.\tPASS\tEND=2300;SVTYPE=DEL;SVLEN=300;AF=0.4;ALGORITHMS=manta;EVIDENCE=PE,SR;PREDICTED_LOF=PRKN\tGT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ\t0/1:2:PE,SR:90:90:1:2",
                # Depth-only but strong RD_GQ/OGQ — kept
                "chr6\t3000\tstrong_depth\tN\t<DEL>\t.\tPASS\tEND=3300;SVTYPE=DEL;SVLEN=300;AF=0.4;ALGORITHMS=depth;EVIDENCE=RD;PREDICTED_LOF=PRKN\tGT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ\t0/1:2:RD:80:20:1:20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    panel = pd.DataFrame({"gene": ["SOX6", "PRKN", "OPA1"], "gene_id": ["1", "2", "3"]})
    table, _ = sv.ingest_sv(vcf, genes, exons, panel, snv_samples=["SAMPLE_A"], auto_index=False)
    ids = set(table["sv_id"].astype(str))
    assert "sox6_art" in ids
    assert "opa1_art" in ids
    assert "pesr" in ids
    assert "strong_depth" in ids
    sox6 = table[table["sv_id"] == "sox6_art"].iloc[0]
    assert bool(sox6["sv_artifact"]) is True
    assert bool(sox6["scored_dosage"]) is False
    kept = table[table["sv_id"] == "strong_depth"].iloc[0]
    assert kept["sv_algorithms"] == "depth"
    assert kept["sv_evidence"] == "RD"
    assert float(kept["rd_gq"]) == 20.0


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


def test_snca_tiny_intronic_dup_is_not_scored(tmp_path):
    genes, exons = _beds(tmp_path)
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Type>",
                "##INFO=<ID=PREDICTED_INTRONIC,Number=.,Type=String,Description=Intron>",
                "##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description=SR>",
                "##INFO=<ID=FREQ_HET,Number=1,Type=Float,Description=Freq>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=CN>",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=ECN>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_C",
                "chr4\t1500\ttiny\tN\t<DUP>\t.\tPASS\tEND=1624;SVTYPE=DUP;PREDICTED_INTRONIC=SNCA;HIGH_SR_BACKGROUND;FREQ_HET=0.1\tGT:CN:ECN\t0/1:3:2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    panel = pd.DataFrame({"gene": ["SNCA"], "gene_id": ["1"]})
    table, _ = sv.ingest_sv(vcf, genes, exons, panel, snv_samples=["SAMPLE_C"], auto_index=False)
    row = table.iloc[0]
    assert row["gene"] == "SNCA"
    assert bool(row["scored_dosage"]) is False
    assert row["dosage_mechanism"] == "other"
    assert row["consequence"] != "copy_number_gain"
    assert bool(row["high_sr_background"]) is True
    assert float(row["freq_het"]) == pytest.approx(0.1)


def test_snca_tiny_copy_gain_prediction_is_not_scored(tmp_path):
    genes, exons = _beds(tmp_path)
    vcf = tmp_path / "sv.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=Type>",
                "##INFO=<ID=PREDICTED_COPY_GAIN,Number=.,Type=String,Description=Gain>",
                "##INFO=<ID=ClinVar_CLNSIG,Number=1,Type=String,Description=ClinVar>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=CN>",
                "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=ECN>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_C",
                "chr4\t3000\ttinygain\tN\t<DUP>\t.\tPASS\tEND=3124;SVTYPE=DUP;PREDICTED_COPY_GAIN=SNCA\tGT:CN:ECN\t0/1:3:2",
                "chr4\t900\tplp\tN\t<DUP>\t.\tPASS\tEND=20000;SVTYPE=DUP;PREDICTED_COPY_GAIN=SNCA;ClinVar_CLNSIG=Pathogenic\tGT:CN:ECN\t0/1:4:2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    panel = pd.DataFrame({"gene": ["SNCA"], "gene_id": ["1"]})
    table, _ = sv.ingest_sv(vcf, genes, exons, panel, snv_samples=["SAMPLE_C"], auto_index=False)
    tiny = table[table["sv_id"] == "tinygain"].iloc[0]
    assert bool(tiny["scored_dosage"]) is False
    assert tiny["consequence"] != "copy_number_gain"
    scored = table[table["sv_id"] == "plp"].iloc[0]
    assert bool(scored["scored_dosage"]) is True
    assert scored["consequence"] == "copy_number_gain"
    assert scored["clinvar_class"] == "P/LP"


def test_str_expansion_sets_classification():
    ad = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "FGF14",
                "str_locus": "FGF14",
                "inheritance_status": "str_expanded",
                "str_expanded": True,
                "str_n_expanded_alleles": 1,
                "str_motif_pathogenic": True,
                "variant_id": "fgf14",
            }
        ]
    )
    out = liability.classify_cohort(ad, ["S1"], cfg={"str_ar_loci": ["RFC1", "CSTB", "FXN"]})
    assert out.iloc[0]["classification"] == "solved"

    rfc1 = pd.DataFrame(
        [
            {
                "sample_id": "S2",
                "gene": "RFC1",
                "str_locus": "RFC1",
                "inheritance_status": "str_expanded",
                "str_expanded": True,
                "str_n_expanded_alleles": 1,
                "str_motif_pathogenic": False,
                "variant_id": "rfc1",
            }
        ]
    )
    out = liability.classify_cohort(rfc1, ["S2"], cfg={"str_ar_loci": ["RFC1", "CSTB", "FXN"]})
    assert out.iloc[0]["classification"] == "partially_explained"

    biallelic = rfc1.copy()
    biallelic["sample_id"] = "S3"
    biallelic["str_n_expanded_alleles"] = 2
    biallelic["str_motif_pathogenic"] = True
    out = liability.classify_cohort(biallelic, ["S3"], cfg={"str_ar_loci": ["RFC1", "CSTB", "FXN"]})
    assert out.iloc[0]["classification"] == "solved"


def test_intron_only_ad_is_partial_including_mapt():
    intron = pd.DataFrame(
        [
            {
                "sample_id": "S1",
                "gene": "GRN",
                "inheritance_status": "AD_candidate",
                "consequence": "intron_variant",
                "clinvar_class": "none",
                "mechanism": "lof",
                "variant_id": "grn-intron",
            }
        ]
    )
    reduced = {"reduced_penetrance_genes": ["GBA1", "LRRK2"]}
    out = liability.classify_cohort(intron, ["S1"], cfg=reduced)
    assert out.iloc[0]["classification"] == "partially_explained"

    mapt = intron.copy()
    mapt["gene"] = "MAPT"
    mapt["mechanism"] = "splice_gof"
    mapt["consequence"] = "splice_region_variant"
    out = liability.classify_cohort(mapt, ["S1"], cfg=reduced)
    assert out.iloc[0]["classification"] == "partially_explained"

    donor = mapt.copy()
    donor["consequence"] = "splice_donor_variant"
    out = liability.classify_cohort(donor, ["S1"], cfg=reduced)
    assert out.iloc[0]["classification"] == "solved"

    plp = mapt.copy()
    plp["clinvar_class"] = "P/LP"
    out = liability.classify_cohort(plp, ["S1"], cfg=reduced)
    assert out.iloc[0]["classification"] == "solved"


def test_large_vcf_is_not_scanned(tmp_path, monkeypatch):
    from pdvar import vcf_stream

    vcf = tmp_path / "big.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\nchr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n", encoding="utf-8")
    monkeypatch.setattr(vcf_stream, "STREAM_MAX_BYTES", 1)
    with pytest.raises(RuntimeError, match="Refusing to scan"):
        list(vcf_stream.iter_vcf(vcf, allow_full_scan=False))
    records = list(vcf_stream.iter_vcf(vcf, allow_full_scan=True))
    assert records[0]["pos"] == 1
