"""Static checks that the entry workflow is the per-sample design."""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Public hg38 chrM bundle from GATK 4.7.0.0 ExampleInputsMitochondriaPipeline.json.
# blacklisted_sites is the BED that tag ships; the public prefix has no VCF.
_MT = "gs://gcp-public-data--broad-references/hg38/v0/chrM"
PUBLIC_MT = {
    "PdSample.mt_dict": f"{_MT}/Homo_sapiens_assembly38.chrM.dict",
    "PdSample.mt_fasta": f"{_MT}/Homo_sapiens_assembly38.chrM.fasta",
    "PdSample.mt_fasta_index": f"{_MT}/Homo_sapiens_assembly38.chrM.fasta.fai",
    "PdSample.mt_amb": f"{_MT}/Homo_sapiens_assembly38.chrM.fasta.amb",
    "PdSample.mt_ann": f"{_MT}/Homo_sapiens_assembly38.chrM.fasta.ann",
    "PdSample.mt_bwt": f"{_MT}/Homo_sapiens_assembly38.chrM.fasta.bwt",
    "PdSample.mt_pac": f"{_MT}/Homo_sapiens_assembly38.chrM.fasta.pac",
    "PdSample.mt_sa": f"{_MT}/Homo_sapiens_assembly38.chrM.fasta.sa",
    "PdSample.blacklisted_sites": f"{_MT}/blacklist_sites.hg38.chrM.bed",
    "PdSample.blacklisted_sites_index": f"{_MT}/blacklist_sites.hg38.chrM.bed.idx",
    "PdSample.mt_shifted_dict": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict",
    "PdSample.mt_shifted_fasta": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta",
    "PdSample.mt_shifted_fasta_index": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai",
    "PdSample.mt_shifted_amb": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb",
    "PdSample.mt_shifted_ann": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann",
    "PdSample.mt_shifted_bwt": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt",
    "PdSample.mt_shifted_pac": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac",
    "PdSample.mt_shifted_sa": f"{_MT}/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa",
    "PdSample.shift_back_chain": f"{_MT}/ShiftBack.chain",
    "PdSample.control_region_shifted_reference_interval_list": f"{_MT}/control_region_shifted.chrM.interval_list",
    "PdSample.non_control_region_interval_list": f"{_MT}/non_control_region.chrM.interval_list",
}

# Documented in workflow/terra/STAGING.md. The only allowed WORKSPACE_BUCKET form.
STAGED = {
    "PdSample.vep_cache": "gs://YOUR_WORKSPACE_BUCKET/resources/vep/homo_sapiens_vep_116_GRCh38.tar.gz",
    "PdSample.somalier_sites": "gs://YOUR_WORKSPACE_BUCKET/resources/somalier/sites.hg38.vcf.gz",
}
EXAMPLE_SAMPLE = {
    "PdSample.cram": "gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.cram",
    "PdSample.crai": "gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.cram.crai",
}


def test_pdsample_is_the_only_entry_workflow():
    text = (ROOT / "workflow" / "wdl" / "PdSample.wdl").read_text(encoding="utf-8")
    assert "workflow PdSample" in text
    for banned in (
        "force_rerun",
        "dry_run",
        "PdCohort",
        "JointGenotyping",
        "results_prefix",
        "submissions_prefix",
        "allow_fallback",
        "sep(",
    ):
        assert banned not in text
    assert "gatk-sv/v0.29-beta/wdl/GATKSVPipelineSingleSample.wdl" in text
    assert "gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl" in text
    dockstore = (ROOT / ".dockstore.yml").read_text(encoding="utf-8")
    assert "PdSample.wdl" in dockstore
    assert "PdCohort" not in dockstore
    assert "gatk_sv_ready" not in text
    assert "WriteSvStatus" not in text
    example = (ROOT / "workflow" / "terra" / "PdSample.example.json").read_text(encoding="utf-8")
    assert "SYNTHETIC_SAMPLE" in example
    assert "gs://fc-" not in example
    assert "raw.githubusercontent.com/gqworks/pdsample-terra/main/" in example
    assert "ref_panel_1kg_v2-contig-ploidy-model.tar.gz" in example
    terra = (ROOT / "workflow" / "terra" / "PdSample.terra_inputs.json").read_text(encoding="utf-8")
    assert '"PdSample.sample_id": "${this.sample_id}"' in terra
    assert '"PdSample.cram": "${this.cram}"' in terra
    assert '"PdSample.crai": "${this.crai}"' in terra
    assert '"PdSample.sex": "${this.sex}"' in terra
    assert '"PdSample.sv_ref_panel_vcf": "${workspace.ref_panel_vcf}"' in terra
    assert '"PdSample.sv_empty_file": "${workspace.reference_empty_file}"' in terra
    assert "gs://fc-" not in terra
    assert "workspace.sv_pipeline_docker" not in terra
    assert "workspace.manta_docker" not in terra


def test_json_placeholders_are_public_or_documented():
    staging = (ROOT / "workflow" / "terra" / "STAGING.md").read_text(encoding="utf-8")
    assert "gs://fc-" not in staging
    for path in (
        ROOT / "workflow" / "terra" / "PdSample.example.json",
        ROOT / "workflow" / "terra" / "PdSample.terra_inputs.json",
    ):
        text = path.read_text(encoding="utf-8")
        data = json.loads(text)
        assert "gs://fc-" not in text
        assert "REPLACE_" not in text
        assert re.search(r"(?<!YOUR_)WORKSPACE_BUCKET", text) is None
        for key, value in PUBLIC_MT.items():
            assert data[key] == value
        for key, value in STAGED.items():
            assert data[key] == value
            assert value in staging
            assert key in staging
        for value in data.values():
            if isinstance(value, str) and "YOUR_WORKSPACE_BUCKET" in value:
                assert value in staging
    example = json.loads(
        (ROOT / "workflow" / "terra" / "PdSample.example.json").read_text(encoding="utf-8")
    )
    terra = json.loads(
        (ROOT / "workflow" / "terra" / "PdSample.terra_inputs.json").read_text(encoding="utf-8")
    )
    for key, value in EXAMPLE_SAMPLE.items():
        assert example[key] == value
        assert value in staging
        assert key in staging
    assert terra["PdSample.cram"] == "${this.cram}"
    assert terra["PdSample.crai"] == "${this.crai}"
    assert "YOUR_WORKSPACE_BUCKET" not in terra["PdSample.cram"]
