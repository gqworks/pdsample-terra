"""MitImpact CSV join onto an mtDNA VCF."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "annotate_mtdna.py"


def _load_mod():
    import importlib.util

    spec = importlib.util.spec_from_file_location("annotate_mtdna", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_annotate_mtdna_passthrough_without_table(tmp_path):
    annotate = _load_mod()
    vcf = tmp_path / "mt.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "chrM\t11778\t.\tG\tA\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    status = tmp_path / "status.txt"
    annotate.annotate(vcf, out, None, status)
    assert status.read_text(encoding="utf-8").strip() == "not_assessed"
    assert "chrM\t11778" in out.read_text(encoding="utf-8")
    assert "MI_APOGEE2" not in out.read_text(encoding="utf-8")


def test_annotate_mtdna_joins_mitimpact_stub(tmp_path):
    annotate = _load_mod()
    table = tmp_path / "mitimpact.csv"
    table.write_text(
        "Start,Ref,Alt,Gene_symbol,class,ref_aa,alt_aa,aa_pos,HGVS_genomic,"
        "APOGEE2,APOGEE2_score,Clinvar_CLNSIG,Clinvar_CLNDN,MITOMAP_Disease_Status,"
        "MITOMAP_Disease_GenBank_Freq,MLC_score,"
        "AF_hom_gnomad,AF_het_gnomad,MitoTIP_Score,Mitotip_Quartile\n"
        "11778,G,A,MT-ND4,protein-coding,R,H,340,NC_012920.1:m.11778G>A,"
        "Pathogenic,0.926,Pathogenic,Leber_optic_atrophy,Cfrm-[P],"
        "0.00012,0.72,"
        "0.0002,0.0002,,\n",
        encoding="utf-8",
    )
    vcf = tmp_path / "mt.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "chrM\t11778\t.\tG\tA\t.\tPASS\tGENE=MT-ND4\tGT\t0/1",
                "chrM\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    status = tmp_path / "status.txt"
    annotate.annotate(vcf, out, table, status)
    text = out.read_text(encoding="utf-8")
    assert status.read_text(encoding="utf-8").strip() == "called"
    assert "MI_APOGEE2=Pathogenic" in text
    assert "MI_ClinVar_CLNSIG=Pathogenic" in text
    assert "MI_MITOMAP_status=Cfrm-[P]" in text
    assert "MI_gnomAD_AF_hom=0.0002" in text
    assert "MI_MLC_score=0.72" in text
    assert "MI_MITOMAP_GB_freq=0.00012" in text
    assert "MI_gene=MT-ND4" in text
    # Unmatched allele keeps original INFO.
    assert "chrM\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1" in text
