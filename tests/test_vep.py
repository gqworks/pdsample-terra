"""VEP INFO parsing for gnomAD, ClinVar, CADD/REVEL, and in-silico predictors."""

from pathlib import Path

from pdvar.vep import _clinvar_from_info, _gnomad_from_info, _in_silico_from_csq, parse_vep


def test_gnomad_from_info_prefers_max_across_genomes_and_exomes():
    af, pop = _gnomad_from_info(
        {
            "gnomADg_AF": "0.001",
            "gnomADe_AF": "0.02",
            "gnomADg_AF_popmax": "0.003",
            "gnomADe_AF_popmax": "0.01",
        }
    )
    assert af == 0.02
    assert pop == 0.01


def test_gnomad_from_info_uses_lifted_unified_keys():
    af, pop = _gnomad_from_info({"gnomAD_AF": "0.05", "AF_popmax": "0.07"})
    assert af == 0.05
    assert pop == 0.07


def test_gnomad_from_info_legacy_pipe():
    af, pop = _gnomad_from_info({"gnomAD": "0.004|0.006|1"})
    assert af == 0.004
    assert pop == 0.006


def test_clinvar_from_info_reads_clnsig_field():
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
            "FATHMM_score": "-2.1",
            "FATHMM_pred": "D",
        }
    )
    assert dbnsfp["sift"] == 0.02
    assert dbnsfp["sift_pred"] == "d"
    assert dbnsfp["polyphen"] == 0.9
    assert dbnsfp["eigen"] == 12.3
    assert dbnsfp["muttaster"] == 0.99
    assert dbnsfp["muttaster_pred"] == "D"
    assert dbnsfp["fathmm"] == -2.1
    assert dbnsfp["fathmm_pred"] == "D"


def test_parse_vep_fills_cadd_revel_clinvar_class_and_gnomad(tmp_path: Path):
    csq = (
        "Allele|Consequence|SYMBOL|HGVSc|HGVSp|CANONICAL|CADD_PHRED|REVEL|"
        "SIFT|PolyPhen|Eigen-phred_coding|MutationTaster_pred|FATHMM_pred"
    )
    path = tmp_path / "panel.vep.vcf"
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {csq}">',
                '##INFO=<ID=gnomADg_AF,Number=A,Type=Float,Description="gnomAD genomes AF">',
                '##INFO=<ID=gnomADe_AF,Number=A,Type=Float,Description="gnomAD exomes AF">',
                '##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="ClinVar CLNSIG">',
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "1\t100\t.\tA\tG\t.\tPASS\t"
                "CSQ=G|missense_variant|PINK1|c.1A>G|p.X|YES|28.5|0.85|"
                "deleterious(0.01)|probably_damaging(0.99)|11.5|D|D;"
                "gnomADg_AF=0.001;gnomADe_AF=0.01;ClinVar_CLNSIG=Pathogenic\t"
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
    assert row["clinvar"] == "Pathogenic"
    assert row["clinvar_class"] == "P/LP"
    assert row["sift_pred"] == "deleterious"
    assert row["polyphen_pred"] == "probably_damaging"
    assert row["eigen"] == 11.5
    assert row["muttaster_pred"] == "D"
    assert row["fathmm_pred"] == "D"
