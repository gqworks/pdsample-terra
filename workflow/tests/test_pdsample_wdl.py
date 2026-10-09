"""Static checks that the entry workflow is the per-sample design."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


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
