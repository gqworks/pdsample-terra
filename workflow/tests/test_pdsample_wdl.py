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
    "PdSample.introme_gtf": "gs://YOUR_WORKSPACE_BUCKET/resources/gencode/gencode.v44.annotation.gtf.gz",
}
EXAMPLE_SAMPLE = {
    "PdSample.bam": "gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.bam",
    "PdSample.bai": "gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.bam.bai",
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
    assert 'import "vendor/gatk-sv/GATKSVPipelineSingleSample.wdl"' in text
    assert "gatk-sv/v1.1.1/wdl/" in (
        ROOT / "workflow" / "wdl" / "vendor" / "gatk-sv" / "GATKSVPipelineSingleSample.wdl"
    ).read_text(encoding="utf-8")
    assert "gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl" in text
    sample_tasks = (ROOT / "workflow" / "wdl" / "tasks" / "SampleTasks.wdl").read_text(encoding="utf-8")
    assert "pysam==0.24.1" in sample_tasks
    assert "alphagenome==" not in sample_tasks
    assert "ALPHAGENOME_API_KEY" not in sample_tasks
    assert "File? avi_table" in text
    assert "File? avi_table_index" in text
    assert "alphagenome_api_key" not in text
    assert "AnnotateAvi" in text
    dockstore = (ROOT / ".dockstore.yml").read_text(encoding="utf-8")
    assert "PdSample.wdl" in dockstore
    assert "PdCohort" not in dockstore
    assert "gatk_sv_ready" not in text
    assert "WriteSvStatus" not in text
    example = (ROOT / "workflow" / "terra" / "PdSample.example.json").read_text(encoding="utf-8")
    assert "SYNTHETIC_SAMPLE" in example
    assert "gs://fc-" not in example
    assert "raw.githubusercontent.com/gqworks/pdsample-terra/main/" not in example
    assert "gatk-sv-ref-panel-1kg-v1-0-1" in example
    assert "PdSample.sv_ref_std_scramble_vcf_tar" in example
    assert "PdSample.sv_empty_file" not in example
    assert '"../' not in example
    assert "alphagenome_api_key" not in example
    assert "PdSample.avi_table" not in example
    terra = (ROOT / "workflow" / "terra" / "PdSample.terra_inputs.json").read_text(encoding="utf-8")
    assert '"PdSample.sample_id": "${this.sample_id}"' in terra
    assert '"PdSample.bam": "${this.bam_or_cram_file}"' in terra
    assert '"PdSample.bai": "${this.bai}"' in terra
    assert '"PdSample.sex": "${this.sex}"' in terra
    assert '"PdSample.age_onset": "${this.age_onset}"' in terra
    assert '"PdSample.symptoms": "${this.symptoms}"' in terra
    assert '"PdSample.ethnicity": "${this.ethnicity}"' in terra
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "pheno[Phenotype report or age, symptoms, ancestry] --> assess" in readme
    assert '"PdSample.phenotype_report": "${this.phenotype_report}"' in terra
    assert "File? phenotype_report" in text
    tasks = (ROOT / "workflow" / "wdl" / "tasks" / "SampleTasks.wdl").read_text(encoding="utf-8")
    assert "PHENOTYPE_REPORT" in tasks
    assert "String? ethnicity" in text
    assert "ethnicity = ethnicity" in text
    assert "write_json" in tasks
    assert '"ethnicity": ethnicity_value' in tasks
    assert '"PdSample.sv_ref_panel_vcf": "${workspace.ref_panel_vcf}"' in terra
    assert '"PdSample.sv_ref_std_scramble_vcf_tar": "${workspace.ref_panel_std_scramble_vcf_tar}"' in terra
    assert "sv_empty_file" not in terra
    assert "gs://fc-" not in terra
    assert "workspace.sv_pipeline_docker" not in terra
    assert "workspace.manta_docker" not in terra
    assert "raw.githubusercontent.com" not in terra
    assert "check_reference_py" not in terra


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
    assert terra["PdSample.bam"] == "${this.bam_or_cram_file}"
    assert terra["PdSample.bai"] == "${this.bai}"
    assert "YOUR_WORKSPACE_BUCKET" not in terra["PdSample.bam"]


def test_wdl_file_defaults_are_not_relative_paths():
    """Terra localizes every File, including defaults. Relative paths fail there."""
    wdl_root = ROOT / "workflow" / "wdl"
    for path in sorted(wdl_root.rglob("*.wdl")):
        if "stubs" in path.parts:
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.startswith("File ") and ("\"../" in stripped or "\"../../" in stripped):
                raise AssertionError(f"{path}: relative File default: {stripped}")


def test_example_json_has_no_relative_file_paths():
    """Terra GCP Batch only localizes gs://; relative paths in JSON are passed literally."""
    example = (ROOT / "workflow" / "terra" / "PdSample.example.json").read_text(encoding="utf-8")
    assert '"../' not in example
    assert '"../../' not in example


def test_file_inputs_do_not_use_https_on_terra_json():
    """Terra GCP Batch localizes gs:// only; https File paths fail at runtime."""
    fileish_prefixes = (
        "PdSample.panel_",
        "PdSample.gene_",
        "PdSample.exon_",
        "PdSample.eh_",
        "PdSample.config_",
        "PdSample.columns_",
        "PdSample.phenotypes_",
        "PdSample.check_",
        "PdSample.qc_",
        "PdSample.annotate_",
        "PdSample.verifybamid_",
        "PdSample.introme_gtf",
        "PdSample.vep_cache",
        "PdSample.somalier_",
    )
    terra = json.loads(
        (ROOT / "workflow" / "terra" / "PdSample.terra_inputs.json").read_text(encoding="utf-8")
    )
    for key, value in terra.items():
        if not any(key.startswith(p) for p in fileish_prefixes):
            continue
        if isinstance(value, str) and value.startswith("https://"):
            raise AssertionError(f"{key} must be gs:// on Terra, not {value}")


def test_cp_uses_localized_file_interpolation():
    """Coercing File to String in WDL leaves gs:// paths; shell cp cannot use them."""
    wdl_root = ROOT / "workflow" / "wdl"
    banned = (
        "dict_path",
        'cp "~{supplied_tar}"',
    )
    for path in sorted(wdl_root.rglob("*.wdl")):
        if "stubs" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        for fragment in banned:
            assert fragment not in text, f"{path}: found banned fragment {fragment!r}"
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("cp ") and "~{" in stripped:
                assert "_path}" not in stripped, f"{path}: cp uses string-coerced *_path: {stripped}"
        for line in text.splitlines():
            if '"" + select_first' in line:
                assert False, f"{path}: File coerced to String: {line.strip()}"


def test_alignment_readers_accept_picard_bai_names():
    """htslib looks for sample.bam.bai; Picard writes sample.bai."""
    text = (ROOT / "workflow" / "wdl" / "tasks" / "SampleTasks.wdl").read_text(encoding="utf-8")
    assert text.count('ln -s "$index" "${align}.bai"') == 5


def test_requester_pays_flag_is_a_haplotypecaller_argument():
    """GATK treats a token before the tool name as the command."""
    text = (ROOT / "workflow" / "wdl" / "tasks" / "SampleTasks.wdl").read_text(encoding="utf-8")
    assert "HaplotypeCaller \\\n      $rp_arg \\\n" in text


def test_mito_and_hc_prefer_original_ref_dict_for_adjacency():
    """Imported mito (and HC -R) need <fasta>.dict co-localized; renamed reference.dict breaks that."""
    text = (ROOT / "workflow" / "wdl" / "PdSample.wdl").read_text(encoding="utf-8")
    assert "File gatk_ref_dict = select_first([ref_dict, select_first(reference_dict_files)])" in text
    assert text.count("ref_dict = gatk_ref_dict") == 2
    # GenotypeHardFilter keeps the copied dict and bundles beside -R in-task.
    assert "call Tasks.GenotypeHardFilter" in text
    ghf = text.split("call Tasks.GenotypeHardFilter", 1)[1].split("call Tasks.NormPass", 1)[0]
    assert "ref_dict = select_first(reference_dict_files)," in ghf


def test_genotype_hard_filter_bundles_reference_dict_with_fasta():
    """GenotypeGVCFs looks for <ref>.dict next to the FASTA path passed to -R."""
    text = (ROOT / "workflow" / "wdl" / "tasks" / "SampleTasks.wdl").read_text(encoding="utf-8")
    assert text.count('dict_name="${ref_base%.*}.dict"') == 1
    assert 'GenotypeGVCFs -R "$ref_fasta_path"' in text
    assert "HaplotypeCaller \\\n      $rp_arg \\\n      -R ~{ref_fasta}" in text
    assert "readlink -f ~{ref_fasta}" not in text


def test_qc_tasks_size_disk_for_the_full_alignment():
    """Tasks that copy the BAM must budget the whole file, its index, and headroom."""
    text = (ROOT / "workflow" / "wdl" / "tasks" / "SampleTasks.wdl").read_text(encoding="utf-8")
    assert 'size(bam, "GiB") * 0.2' not in text
    assert text.count('size(bam, "GiB") * 2 + size(bam_index, "GiB")') >= 5


def test_sv_alignment_is_released_as_bam_after_qc():
    """Manta 1.6 needs a .bam filename so it loads sample.bam.bai. BAM-only inputs."""
    workflow = (ROOT / "workflow" / "wdl" / "PdSample.wdl").read_text(encoding="utf-8")
    tasks = (ROOT / "workflow" / "wdl" / "tasks" / "SampleTasks.wdl").read_text(encoding="utf-8")
    assert "task StageIndexedAlignment" not in tasks
    assert "StageIndexedAlignment" not in workflow
    assert "call Tasks.ReleaseAfterQc as ReleaseBam" in workflow
    assert "call Tasks.ReleaseAfterQc as ReleaseBai" in workflow
    assert 'name = sample_id + ".bam"' in workflow
    assert 'name = sample_id + ".bam.bai"' in workflow
    assert "bam_or_cram_file = ReleaseBam.released" in workflow
    assert "bam_or_cram_index = ReleaseBai.released" in workflow
    assert "File bam" in workflow
    assert "File bai" in workflow
    assert "File cram" not in workflow
    assert "task ExtractBamHeader" in tasks
    assert "*.cram)" not in tasks


def test_mito_leaves_initial_filter_vcf_uncompressed():
    """SplitMultiAllelics reads the InitialFilter VCF without its index."""
    text = (ROOT / "workflow" / "wdl" / "PdSample.wdl").read_text(encoding="utf-8")
    assert "compress_output_vcf = true" not in text
    assert "compress_output_vcf = false" in text


def test_introme_indexes_gtf_as_gff_and_guards_missing_strand():
    """ESE scoring reads INFO/strand. .gtf.gz is not a tabix preset name."""
    text = (ROOT / "workflow" / "wdl" / "tasks" / "IntromeTasks.wdl").read_text(encoding="utf-8")
    assert "tabix -f sorted.annotation.gtf.gz" not in text
    assert "tabix -f -p gff sorted.annotation.gtf.gz" in text
    assert "tabix -f -p gff \"$gtf_path\"" in text
    assert '##INFO=<ID=strand,Number=.,Type=String,Description="GENCODE gene strand">' in text
    assert "raw_strand = record.info['strand']" in text
    assert "if raw_strand is None:" in text


def test_introme_docker_tags_are_published():
    """Unpublished tags fail the Batch pull before the command, so Terra has no log."""
    introme = (ROOT / "workflow" / "wdl" / "Introme.wdl").read_text(encoding="utf-8")
    sample = (ROOT / "workflow" / "wdl" / "PdSample.wdl").read_text(encoding="utf-8")
    tasks = (ROOT / "workflow" / "wdl" / "tasks" / "IntromeTasks.wdl").read_text(encoding="utf-8")
    for text in (introme, sample, tasks):
        assert "mmsplice:v2.2" not in text
        assert "spliceogen:v2.1" not in text
    assert 'mmsplice_cpu_docker = "headoncollusion/mmsplice:v2.1-cpu"' in introme
    assert 'mmsplice_gpu_docker = "headoncollusion/mmsplice:v2.1-gpu"' in introme
    assert 'spliceogen_docker = "headoncollusion/spliceogen:v2.0"' in introme
    assert 'mmsplice_cpu_docker = "headoncollusion/mmsplice:v2.1-cpu"' in sample
    assert 'spliceogen_docker = "headoncollusion/spliceogen:v2.0"' in sample


def test_introme_gpu_attributes_are_absent_unless_a_gpu_is_requested():
    """gpuCount: 0 fails in the GCP Batch job actor, before any stderr exists."""
    text = (ROOT / "workflow" / "wdl" / "tasks" / "IntromeTasks.wdl").read_text(encoding="utf-8")
    assert "gpuCount: gpu_count" not in text.split("gcp:", 1)[0]
    for name, nxt in (
        ("SpliceAI", "MMSplice"),
        ("MMSplice", "Pangolin"),
        ("Pangolin", "Spip"),
    ):
        block = text.split(f"task {name} {{", 1)[1].split(f"task {nxt} {{", 1)[0]
        runtime = block.split("runtime {", 1)[1]
        assert "\n    gpuCount:" not in runtime
        assert "\n    gpuType:" not in runtime
        assert "\n    nvidiaDriverVersion:" not in runtime
        assert "if (use_gpu && gpu_count > 0) then object {" in runtime
        assert "gpuCount: gpu_count," in runtime
        assert "gpuType: gpu_type," in runtime
        assert "nvidiaDriverVersion: nvidia_driver_version" in runtime
        assert "} else object {}" in runtime


def test_introme_cpu_tools_are_given_a_thread_count():
    """Neural nets have no CLI thread flag. They read OMP and TensorFlow env vars."""
    tasks = (ROOT / "workflow" / "wdl" / "tasks" / "IntromeTasks.wdl").read_text(encoding="utf-8")
    workflow = (ROOT / "workflow" / "wdl" / "Introme.wdl").read_text(encoding="utf-8")
    assert tasks.count('export TF_NUM_INTRAOP_THREADS="$threads"') == 2
    assert tasks.count('export OMP_NUM_THREADS="$threads"') >= 4
    assert 'exec /usr/bin/sort --parallel=${threads}' in tasks
    assert "Int pangolin_cpu = if use_gpu then 1 else 16" in workflow
    assert "Int spliceai_cpu = if use_gpu then 1 else 16" in workflow
    assert "Int mmsplice_cpu = if use_gpu then 1 else 16" in workflow
    assert "memory_gb = spliceai_mem" in workflow
    assert "memory_gb = mmsplice_mem" in workflow


def test_spliceogen_scratch_is_on_the_task_disk():
    """RUN.sh writes under /Spliceogen, which is the boot disk unless redirected."""
    text = (ROOT / "workflow" / "wdl" / "tasks" / "IntromeTasks.wdl").read_text(encoding="utf-8")
    spliceogen = text.split("task Spliceogen {", 1)[1].split("task IntromeFunctions {", 1)[0]
    assert 'ln -s "$PWD/scratch/data" /Spliceogen/data' in spliceogen
    assert 'ln -s "$PWD/scratch/temp" /Spliceogen/temp' in spliceogen
    assert 'ln -s "$PWD/scratch/output" /Spliceogen/output' in spliceogen
    assert 'export TMPDIR="$PWD/scratch/temp"' in spliceogen
    assert "Int preemptible = 0" in spliceogen
    assert "Int memory_gb = 16" in spliceogen
    assert 'export JAVA_TOOL_OPTIONS="-Xmx8g"' in spliceogen
    assert 'VCF_PATH="$(realpath ~{vcf})"' in spliceogen
    assert 'GTF_PATH="$(realpath ~{gtf})"' in spliceogen
    assert 'FASTA_PATH="$ORIG_DIR/ref/ref.fa"' in spliceogen
    assert '-input "$VCF_PATH"' in spliceogen
    assert '-fasta "$FASTA_PATH"' in spliceogen
    assert '-gtf "$GTF_PATH"' in spliceogen
    assert 'realpath "$ORIG_DIR/~{vcf}"' not in spliceogen
    assert 'realpath "$ORIG_DIR/~{gtf}"' not in spliceogen


def test_gatksv_callers_use_sixteen_cores():
    """Manta and Wham read cpu_cores. Upstream defaults are lower."""
    text = (ROOT / "workflow" / "wdl" / "PdSample.wdl").read_text(encoding="utf-8")
    assert "Float sv_manta_jobs_per_cpu = 1.0" in text
    assert "Int sv_caller_cpu_cores = 16" in text
    assert "Float sv_caller_mem_gb = 32.0" in text
    assert "manta_jobs_per_cpu = sv_manta_jobs_per_cpu" in text
    assert "runtime_attr_manta = object {" in text
    assert "runtime_attr_wham = object {" in text
    assert "cpu_cores: sv_caller_cpu_cores" in text
    assert "mem_gb: sv_caller_mem_gb" in text


def test_gatksv_uses_v1_1_1_with_scramble_not_melt():
    """Dockstore SingleSamplePipeline v1.1.1: Scramble on, MELT off, new panel inputs."""
    text = (ROOT / "workflow" / "wdl" / "PdSample.wdl").read_text(encoding="utf-8")
    vendor = (
        ROOT / "workflow" / "wdl" / "vendor" / "gatk-sv" / "GATKSVPipelineSingleSample.wdl"
    ).read_text(encoding="utf-8")
    assert "Boolean sv_use_melt = false" in text
    assert "Boolean sv_use_scramble = true" in text
    assert "use_melt = sv_use_melt" in text
    assert "use_scramble = sv_use_scramble" in text
    assert "ref_std_scramble_vcf_tar = sv_ref_std_scramble_vcf_tar" in text
    assert "gq_recalibrator_model_file = sv_gq_recalibrator_model_file" in text
    assert "sl_filter_args = sv_sl_filter_args" in text
    assert "HERVK_reference = sv_hervk_reference" in text
    assert "LINE1_reference = sv_line1_reference" in text
    assert "refine_complex_n_per_split = sv_refine_complex_n_per_split" in text
    assert "empty_file" not in text
    assert "HERVK_reference=HERVK_reference" in vendor
    assert "n_per_split=refine_complex_n_per_split" in vendor


def test_runtime_disk_gb_is_integer_ceiled():
    """GCP Batch rejects fractional local-disk sizes (e.g. 13.026... HDD)."""
    wdl_root = ROOT / "workflow" / "wdl"
    for path in sorted(wdl_root.rglob("*.wdl")):
        if "stubs" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        assert "Float disk_gb" not in text, path
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("disk_gb =") or stripped.startswith("Int disk_gb ="):
                assert "ceil(" in stripped, f"{path}: disk_gb must use ceil(): {stripped}"
                assert stripped.startswith("Int disk_gb"), f"{path}: disk_gb must be Int: {stripped}"
