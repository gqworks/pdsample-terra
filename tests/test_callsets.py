"""Per-sample ExpansionHunter, Gauchian, and mtDNA ingest."""

import json
from pathlib import Path

import pandas as pd

from pdvar import callsets, filters


def _eh(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End>",
                "##INFO=<ID=REPID,Number=1,Type=String,Description=Locus>",
                "##INFO=<ID=REF,Number=1,Type=Integer,Description=Reference repeats>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=REPCN,Number=1,Type=String,Description=Repeat counts>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B",
                "chr4\t39348425\t.\tA\t<STR>\t.\tPASS\tEND=39348479;REPID=RFC1;REF=11\tGT:REPCN\t0/1:11/80\t0/0:11/11",
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
                "##FORMAT=<ID=AF,Number=A,Type=Float,Description=Heteroplasmy>",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=DP>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A",
                "chrM\t11778\t.\tG\tA\t.\tPASS\tGENE=MT-ND4\tGT:AF:DP\t0/1:0.42:800",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_sample_callsets_ingested_and_kept(tmp_path, cfg):
    eh = tmp_path / "eh.vcf"
    gba = tmp_path / "gba.json"
    mt = tmp_path / "mt.vcf"
    _eh(eh)
    gba.write_text(
        '[{"sample_id": "SAMPLE_A", "GBA_copy_number": 0, "recombinant_allele": ""},'
        '{"sample_id": "SAMPLE_B", "GBA_copy_number": 2, "recombinant_allele": "RecNciI", "pathogenic": false}]',
        encoding="utf-8",
    )
    _mt(mt)
    table, assessment = callsets.ingest_sample_callsets(
        "SAMPLE_A", expansionhunter=eh, gauchian=gba, mtdna=mt, cfg=cfg
    )
    assert set(table["hit_source"]) == {"expansionhunter", "gauchian", "mtdna"}
    rfc1 = table[table["str_locus"] == "RFC1"].iloc[0]
    assert rfc1["sample_id"] == "SAMPLE_A"
    assert rfc1["gene"] == "RFC1"
    # Fixture allele 11/80 is below the cited RFC1 AAGGG floor (250).
    assert bool(rfc1["str_expanded"]) is False
    deletion = table[table["hit_source"] == "gauchian"].iloc[0]
    assert deletion["consequence"] == "transcript_ablation"
    assert deletion["zygosity"] == "hom"
    mito = table[table["hit_source"] == "mtdna"].iloc[0]
    assert mito["gene"] == "MT-ND4"
    assert mito["heteroplasmy"] == 0.42
    status = assessment.set_index("module")["callset_status"]
    assert status["expansionhunter"] == "STR carrier"

    other, other_status = callsets.ingest_sample_callsets("SAMPLE_B", expansionhunter=eh, gauchian=gba, cfg=cfg)
    recombinant = other[other["hit_source"] == "gauchian"].iloc[0]
    assert recombinant["consequence"] == "structural_variant"
    assert recombinant["clinvar_class"] == "none"
    by_module = other_status.set_index("module")["callset_status"]
    assert by_module["expansionhunter"] == "no STR"
    assert by_module["mtdna"] == "mtDNA not assessed"

    panel = pd.DataFrame({"gene": ["PRKN"], "risk_allele_exception": [False]})
    kept = filters.apply_filters(table, panel, cfg["filters"])
    assert set(kept["hit_source"]) == {"expansionhunter", "gauchian", "mtdna"}


def test_atxn8os_above_cutoff_is_expanded(tmp_path, cfg):
    path = tmp_path / "eh_atxn8.vcf"
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=End>",
                "##INFO=<ID=REPID,Number=1,Type=String,Description=Locus>",
                "##INFO=<ID=REF,Number=1,Type=Integer,Description=Reference repeats>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "##FORMAT=<ID=REPCN,Number=1,Type=String,Description=Repeat counts>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A",
                "chr13\t70139383\t.\tA\t<STR>\t.\tPASS\tEND=70139428;REPID=ATXN8OS;REF=15\tGT:REPCN\t1/2:9/101",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    table, _ = callsets.ingest_sample_callsets("SAMPLE_A", expansionhunter=path, cfg=cfg)
    row = table[table["str_locus"] == "ATXN8OS"].iloc[0]
    assert row["repeat_count"] == 101
    assert bool(row["str_expanded"]) is True
    assert "NBK1268" in str(row["str_citation"])


def test_absent_callsets_are_not_assessed():
    table, assessment = callsets.ingest_sample_callsets("SAMPLE_A")
    assert table.empty
    assert set(assessment["callset_status"]) == {"STR not assessed", "GBA1 not assessed", "mtDNA not assessed"}


def test_native_gauchian_json_other_variants(tmp_path, cfg):
    """Gauchian 1.0.2 emits {sample: metrics}; other_variants must become GBA1 rows."""
    path = tmp_path / "MH00246.gauchian.json"
    path.write_text(
        json.dumps(
            {
                "MH00246_SYD": {
                    "Coverage_MAD": 0.079,
                    "Median_depth": 38.88,
                    "deletion_CN": 2,
                    "deletion_CN_raw": 1.811,
                    "GBAP1_like_variant_exon9_11": "",
                    "other_variants": "c.681T>G(p.Asn227Lys),c.667T>C(p.Trp223Arg)",
                    "variant_raw_count": {
                        "c.681T>G(p.Asn227Lys)": "11,14",
                        "c.667T>C(p.Trp223Arg)": "10,17",
                    },
                    "haplotypes": {"is_biallelic": False, "is_carrier": False},
                }
            }
        ),
        encoding="utf-8",
    )
    table, assessment = callsets.ingest_sample_callsets("MH00246_SYD", gauchian=path, cfg=cfg)
    assert not table.empty
    assert set(table["hit_source"]) == {"gauchian"}
    assert set(table["gene"]) == {"GBA1"}
    assert set(table["hgvsc"]) == {"c.681T>G", "c.667T>C"}
    assert set(table["clinvar_class"]) == {"P/LP"}
    assert set(table["dosage_event"]) == {"other_variant"}
    c681 = table[table["hgvsc"] == "c.681T>G"].iloc[0]
    assert c681["depth"] == 25.0
    assert abs(float(c681["allele_balance"]) - 11 / 25) < 1e-9
    status = assessment.set_index("module")["callset_status"]
    assert status["gauchian"] == "GBA1 carrier"


def test_drop_redundant_gauchian_snvs_prefers_small_variant():
    df = pd.DataFrame(
        [
            {
                "gene": "GBA1",
                "hit_source": "small_variant",
                "hgvsc": "ENST00000368373.8:c.681T>G",
                "dosage_event": "",
                "variant_id": "1-155238214-A-C",
            },
            {
                "gene": "GBA1",
                "hit_source": "gauchian",
                "hgvsc": "c.681T>G",
                "dosage_event": "other_variant",
                "variant_id": "GBA1-c.681T>G-S",
            },
            {
                "gene": "GBA1",
                "hit_source": "gauchian",
                "hgvsc": "c.667T>C",
                "dosage_event": "other_variant",
                "variant_id": "GBA1-c.667T>C-S",
            },
        ]
    )
    out = callsets.drop_redundant_gauchian_snvs(df)
    assert set(out["variant_id"]) == {"1-155238214-A-C", "GBA1-c.667T>C-S"}


def test_single_foreign_sample_is_relabelled(tmp_path, cfg):
    mt = tmp_path / "mt.vcf"
    _mt(mt)
    table, assessment = callsets.ingest_sample_callsets("SYNTHETIC_SAMPLE", mtdna=mt, cfg=cfg)
    assert table["sample_id"].tolist() == ["SYNTHETIC_SAMPLE"]
    assert assessment.loc[assessment["module"] == "mtdna", "callset_status"].iloc[0] == "mtDNA carrier"
