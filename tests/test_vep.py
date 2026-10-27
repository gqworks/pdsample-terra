"""VEP CSQ parsing for dbNSFP gnomAD/ClinVar/CADD/REVEL and in-silico predictors."""

from pathlib import Path

from pdvar.vep import _clinvar_from_info, _gnomad_from_info, _in_silico_from_csq, parse_vep


def test_gnomad_from_info_prefers_lifted_unified_keys():
    af, pop = _gnomad_from_info({"gnomAD_AF": "0.05", "AF_popmax": "0.07"})
    assert af == 0.05
    assert pop == 0.07


def test_gnomad_from_info_legacy_pipe():
    af, pop = _gnomad_from_info({"gnomAD": "0.004|0.006|1"})
    assert af == 0.004
    assert pop == 0.006


def test_clinvar_from_info_reads_legacy_clnsig_field():
    assert _clinvar_from_info({"ClinVar_CLNSIG": "Pathogenic"}) == "Pathogenic"
    assert _clinvar_from_info({"ClinVar": "Likely_pathogenic|criteria_provided"}) == "Likely_pathogenic"
    assert _clinvar_from_info({}) == ""


def test_in_silico_from_native_and_dbnsfp_fields():
    native = _in_silico_from_csq(
        {"SIFT": "deleterious(0.01)", "PolyPhen": "probably_damaging(0.95)"}
    )
    assert native["sift"] == 0.01
    assert native["sift_pred"] == "deleterious"
    assert native["polyphen"] == 0.95
    assert native["polyphen_pred"] == "probably_damaging"

    dbnsfp = _in_silico_from_csq(
        {
            "SIFT_score": "0.02",
            "SIFT_pred": "D",
            "Polyphen2_HDIV_score": "0.9",
            "Polyphen2_HDIV_pred": "D",
            "Eigen-phred_coding": "12.3",
            "MutationTaster_score": "0.99",
            "MutationTaster_pred": "D",
            "fathmm-XF_coding_score": "0.82",
            "fathmm-XF_coding_pred": "D",
        }
    )
    assert dbnsfp["sift"] == 0.02
    assert dbnsfp["sift_pred"] == "d"
    assert dbnsfp["polyphen"] == 0.9
    assert dbnsfp["eigen"] == 12.3
    assert dbnsfp["muttaster"] == 0.99
    assert dbnsfp["muttaster_pred"] == "D"
    assert dbnsfp["fathmm"] == 0.82
    assert dbnsfp["fathmm_pred"] == "D"

    legacy = _in_silico_from_csq({"FATHMM_score": "-2.1", "FATHMM_pred": "D"})
    assert legacy["fathmm"] == -2.1
    assert legacy["fathmm_pred"] == "D"


def test_parse_vep_reads_dbnsfp_gnomad_clinvar_cadd_revel(tmp_path: Path):
    csq = (
        "Allele|Consequence|SYMBOL|HGVSc|HGVSp|CANONICAL|CADD_phred|REVEL_score|"
        "SIFT|PolyPhen|Eigen-phred_coding|MutationTaster_pred|fathmm-XF_coding_pred|"
        "clinvar_clnsig|gnomAD4.1_joint_AF|gnomAD4.1_joint_POPMAX_AF"
    )
    path = tmp_path / "panel.vep.vcf"
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {csq}">',
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "1\t100\t.\tA\tG\t.\tPASS\t"
                "CSQ=G|missense_variant|PINK1|c.1A>G|p.X|YES|28.5|0.85|"
                "deleterious(0.01)|probably_damaging(0.99)|11.5|D|D|"
                "Pathogenic|0.01|0.02\t"
                "GT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    frame = parse_vep(path, "S1")
    assert len(frame) == 1
    row = frame.iloc[0]
    assert row["cadd"] == 28.5
    assert row["revel"] == 0.85
    assert row["gnomad_af"] == 0.01
    assert row["gnomad_popmax_af"] == 0.02
    assert row["clinvar"] == "Pathogenic"
    assert row["clinvar_class"] == "P/LP"
    assert row["sift_pred"] == "deleterious"
    assert row["polyphen_pred"] == "probably_damaging"
    assert row["eigen"] == 11.5
    assert row["muttaster_pred"] == "D"
    assert row["fathmm_pred"] == "D"
    assert "spliceai" not in frame.columns


def test_parse_vep_accepts_plugin_cadd_revel_names(tmp_path: Path):
    csq = "Allele|Consequence|SYMBOL|HGVSc|HGVSp|CANONICAL|CADD_PHRED|REVEL"
    path = tmp_path / "plugin.vep.vcf"
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {csq}">',
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "1\t200\t.\tC\tT\t.\tPASS\tCSQ=T|missense_variant|PRKN|c.2C>T|p.Y|YES|22.1|0.71\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    row = parse_vep(path, "S1").iloc[0]
    assert row["cadd"] == 22.1
    assert row["revel"] == 0.71
