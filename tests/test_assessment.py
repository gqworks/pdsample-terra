"""SV dosage, expression weights, liability, and biomarkers."""

from pathlib import Path

import pandas as pd
import pytest

from pdvar import biomarkers, isoform, liability, sv
from pdvar.config import load_config
from pdvar.expression_resources import gtex_brain_median
from pdvar.splicing import merge_splicevault, summarize_splice_scores


def _beds(directory: Path) -> tuple[Path, Path]:
    genes = directory / "genes.bed"
    exons = directory / "exons.bed"
    genes.write_text("chr4\t800\t2500\tENSG\tSNCA\nchr6\t900\t4000\tENSG\tPRKN\n", encoding="utf-8")
    exons.write_text("chr4\t1000\t2000\tENSG\tSNCA\nchr6\t1000\t1200\tENSG\tPRKN\n", encoding="utf-8")
    return genes, exons


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
                "chr6\t1000\tdel\tN\t<DEL>\t.\tMULTIALLELIC\tEND=1300;SVTYPE=DEL;PREDICTED_LOF=PRKN\tGT:CN:ECN\t./.:1:2\t0/0:2:2",
                "chr6\t1000\tdrop\tN\t<DEL>\t.\tHIGH_NCR\tEND=1300;SVTYPE=DEL;PREDICTED_LOF=PRKN\tGT:CN:ECN\t0/1:1:2\t0/0:2:2",
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
    prkn = table[(table["gene"] == "PRKN") & (table["sample_id"] == "SAMPLE_A")].iloc[0]
    assert prkn["zygosity"] == "het"
    assert prkn["cn_inferred_gt"] == True  # noqa: E712
    assert prkn["consequence"] == "transcript_ablation"
    assert "drop" not in set(table["sv_id"].astype(str))
    status = assessment.set_index("sample_id")["sv_status"]
    assert status["SAMPLE_D"] == "SV not assessed"
    assert status["SAMPLE_C"] == "SV carrier"


def test_isoform_weight_is_neutral_without_cutoff(tmp_path):
    gct = tmp_path / "tiny.gct"
    gct.write_text(
        "\n".join(
            [
                "#1.2",
                "2\t5",
                "Name\tDescription\tSubstantia nigra\tPutamen\tCaudate\tFrontal Cortex (BA9)\tLiver",
                "ENSG00000185345\tPRKN\t10\t8\t6\t12\t1",
                "ENSG00000145335\tSNCA\t1\t1\t1\t1\t9",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    expression = gtex_brain_median(gct)
    assert expression.loc[expression["gene"] == "PRKN", "gtex_brain_median_tpm"].iloc[0] == 9
    variants = pd.DataFrame(
        [
            {"gene": "PRKN", "consequence": "stop_gained", "chrom": "6", "pos": 10, "evidence_score": 8, "splice_candidate": False},
            {"gene": "SNCA", "consequence": "stop_gained", "chrom": "4", "pos": 10, "evidence_score": 8, "splice_candidate": False},
        ]
    )
    neutral = isoform.weight_lof(variants, expression, {"min_gtex_tpm": None, "unexpressed_weight": 0.25})
    assert set(neutral["lof_weight"]) == {1.0}
    weighted = isoform.weight_lof(variants, expression, {"min_gtex_tpm": 5, "unexpressed_weight": 0.25})
    assert weighted.loc[weighted["gene"] == "SNCA", "lof_weight"].iloc[0] == 0.25
    assert weighted.loc[weighted["gene"] == "PRKN", "lof_weight"].iloc[0] == 1.0
    missing = isoform.weight_lof(variants, None, {})
    assert set(missing["lof_weight_source"]) == {"neutral_missing_resource"}


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
        [{"variant_id": "1-10-A-G", "evidence_score": 4, "evidence_reasons": "spliceai_high", "consequence": "splice_region_variant"}]
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
