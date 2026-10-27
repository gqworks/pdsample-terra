"""VEP INFO parsing for gnomAD genomes/exomes, ClinVar, CADD, and REVEL."""

from pathlib import Path

from pdvar.vep import _clinvar_from_info, _gnomad_from_info, parse_vep


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


def test_parse_vep_fills_cadd_revel_clinvar_class_and_gnomad(tmp_path: Path):
    csq = "Allele|Consequence|SYMBOL|HGVSc|HGVSp|CANONICAL|CADD_PHRED|REVEL"
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
                "CSQ=G|missense_variant|PINK1|c.1A>G|p.X|YES|28.5|0.85;"
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
