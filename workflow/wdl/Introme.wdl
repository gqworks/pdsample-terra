version 1.0

import "tasks/IntromeTasks.wdl" as Tasks

workflow Introme {
  meta {
    description: "Introme v2 infer mode (nextflow_fns d6da48bc). Panel BED and optional rare-variant filter run before the parallel splice tools. No user-built Introme image."
  }
  input {
    File vcf
    File vcf_index
    File ref_fasta
    File ref_fasta_index
    File gtf
    File? bed
    File? introme_assets_tar
    String prefix = "cohort"
    String force_rerun_token = ""
    String genome_build = "hg38"
    String introme_commit = "d6da48bc84abda66f588470f79abcf1244acfb21"
    String archive_url = "https://github.com/CCICB/introme/archive/d6da48bc84abda66f588470f79abcf1244acfb21.tar.gz"

    Boolean restrict_to_panel = true
    Boolean quality_filter = false
    Int min_qual = 200
    Int min_dp = 20
    Int min_ad = 5
    Boolean filter_rare = true
    Float allele_frequency = 0.01
    String af_field = ""

    Boolean use_gpu = false
    String gpu_type = "nvidia-tesla-t4"
    Int gpu_count = 1
    String nvidia_driver_version = "535.183.01"

    String spliceai_db = "gencode.v44.annotation.txt.gz"
    String pangolin_db = "gencode.v44.annotation.db"
    File? spliceai_db_file
    File? pangolin_db_file
    Int spliceai_distance = 1000
    Int spliceai_mask = 0

    String data_preprocessing_docker = "gabyou/data_prepocessing:v2.0"
    String variant_info_docker = "gabyou/variant_info:v3.2"
    String spliceai_cpu_docker = "headoncollusion/spliceai:v1.4-cpu"
    String spliceai_gpu_docker = "headoncollusion/spliceai:v1.4-gpu"
    String mmsplice_cpu_docker = "headoncollusion/mmsplice:v2.2-cpu"
    String mmsplice_gpu_docker = "headoncollusion/mmsplice:v2.2-gpu"
    String pangolin_cpu_docker = "headoncollusion/pangolin:v1.3-cpu"
    String pangolin_gpu_docker = "headoncollusion/pangolin:v1.3-gpu"
    String spip_docker = "headoncollusion/spip:v1.1"
    String spliceogen_docker = "headoncollusion/spliceogen:v2.1"
    String introme_functions_docker = "headoncollusion/ag_check:v1.6"
    String python_docker = "python:3.11.11-slim-bookworm"
    Int preemptible = 2
    String deprecated_introme_docker = ""
  }

  String spliceai_docker = if use_gpu then spliceai_gpu_docker else spliceai_cpu_docker
  String mmsplice_docker = if use_gpu then mmsplice_gpu_docker else mmsplice_cpu_docker
  String pangolin_docker = if use_gpu then pangolin_gpu_docker else pangolin_cpu_docker
  Int gpu_count_eff = if use_gpu then gpu_count else 0
  Int heavy_preemptible = if use_gpu then 0 else 1
  # Upstream CPU profile is 8 CPUs / 2 GB for Pangolin and 3 CPUs / 1 GB for Spliceogen.
  # GCP custom machines need about 0.9 GiB per vCPU, so the CPU path raises those memories.
  Int pangolin_cpu = if use_gpu then 1 else 8
  Int pangolin_mem = if use_gpu then 2 else 8
  Int spliceai_cpu = if use_gpu then 1 else 8
  Int mmsplice_cpu = if use_gpu then 1 else 8

  call Tasks.FetchIntromeAssets {
    input:
      introme_assets_tar = introme_assets_tar,
      introme_commit = introme_commit,
      archive_url = archive_url,
      genome_build = genome_build,
      deprecated_introme_docker = deprecated_introme_docker,
      force_rerun_token = force_rerun_token,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.DataPreprocessing {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gtf = gtf,
      assets_tar = FetchIntromeAssets.assets_tar,
      bed = bed,
      use_bed = restrict_to_panel,
      prefix = prefix,
      force_rerun_token = force_rerun_token,
      docker = data_preprocessing_docker,
      preemptible = preemptible
  }
  if (quality_filter) {
    call Tasks.QualityFilter {
      input:
        vcf = DataPreprocessing.subset_vcf,
        vcf_index = DataPreprocessing.subset_vcf_index,
        prefix = prefix,
        min_qual = min_qual,
        min_dp = min_dp,
        min_ad = min_ad,
        force_rerun_token = force_rerun_token,
        docker = data_preprocessing_docker,
        preemptible = preemptible
    }
  }
  call Tasks.VariantInfo {
    input:
      vcf = select_first([QualityFilter.quality_vcf, DataPreprocessing.subset_vcf]),
      vcf_index = select_first([QualityFilter.quality_vcf_index, DataPreprocessing.subset_vcf_index]),
      gtf = DataPreprocessing.sorted_gtf,
      assets_tar = FetchIntromeAssets.assets_tar,
      prefix = prefix,
      genome_build = genome_build,
      filter_rare = filter_rare,
      allele_frequency = allele_frequency,
      af_field = af_field,
      force_rerun_token = force_rerun_token,
      docker = variant_info_docker,
      preemptible = preemptible
  }
  # Tool stage is scatter-free: Cromwell runs these calls together.
  # Squirls is not here. main.nf leaves that call commented out.
  call Tasks.SpliceAI {
    input:
      vcf = VariantInfo.variant_info_rmanno,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      prefix = prefix,
      genome_build = genome_build,
      spliceai_db = spliceai_db,
      spliceai_db_file = spliceai_db_file,
      distance = spliceai_distance,
      mask = spliceai_mask,
      use_gpu = use_gpu,
      gpu_type = gpu_type,
      gpu_count = gpu_count_eff,
      nvidia_driver_version = nvidia_driver_version,
      force_rerun_token = force_rerun_token,
      docker = spliceai_docker,
      preemptible = heavy_preemptible,
      cpu = spliceai_cpu
  }
  call Tasks.MMSplice {
    input:
      vcf = VariantInfo.variant_info_rmanno,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gtf = DataPreprocessing.sorted_gtf,
      prefix = prefix,
      use_gpu = use_gpu,
      gpu_type = gpu_type,
      gpu_count = gpu_count_eff,
      nvidia_driver_version = nvidia_driver_version,
      force_rerun_token = force_rerun_token,
      docker = mmsplice_docker,
      preemptible = heavy_preemptible,
      cpu = mmsplice_cpu
  }
  call Tasks.Pangolin {
    input:
      vcf = VariantInfo.variant_info_rmanno,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      prefix = prefix,
      genome_build = genome_build,
      pangolin_db = pangolin_db,
      pangolin_db_file = pangolin_db_file,
      use_gpu = use_gpu,
      gpu_type = gpu_type,
      gpu_count = gpu_count_eff,
      nvidia_driver_version = nvidia_driver_version,
      force_rerun_token = force_rerun_token,
      docker = pangolin_docker,
      preemptible = heavy_preemptible,
      cpu = pangolin_cpu,
      memory_gb = pangolin_mem
  }
  call Tasks.Spip {
    input:
      vcf = VariantInfo.variant_info_rmanno,
      prefix = prefix,
      genome_build = genome_build,
      force_rerun_token = force_rerun_token,
      docker = spip_docker,
      preemptible = preemptible
  }
  call Tasks.Spliceogen {
    input:
      vcf = VariantInfo.variant_info_rmanno,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gtf = DataPreprocessing.sorted_gtf,
      prefix = prefix,
      force_rerun_token = force_rerun_token,
      docker = spliceogen_docker,
      preemptible = preemptible
  }
  call Tasks.IntromeFunctions {
    input:
      variant_info_stripped = VariantInfo.variant_info_stripped,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      assets_tar = FetchIntromeAssets.assets_tar,
      prefix = prefix,
      force_rerun_token = force_rerun_token,
      docker = introme_functions_docker,
      preemptible = preemptible
  }
  call Tasks.SplicingAnno {
    input:
      vcf = VariantInfo.variant_info_vcf,
      vcf_index = VariantInfo.variant_info_vcf_index,
      assets_tar = FetchIntromeAssets.assets_tar,
      spliceai_vcf = SpliceAI.spliceai_vcf,
      mmsplice_vcf = MMSplice.mmsplice_vcf,
      pangolin_vcf = Pangolin.pangolin_vcf,
      spip_vcf = Spip.spip_vcf,
      spliceogen_tsv = Spliceogen.spliceogen_tsv,
      ag_check = IntromeFunctions.ag_check,
      ag_check_index = IntromeFunctions.ag_check_index,
      ese_score = IntromeFunctions.ese_score,
      ese_score_index = IntromeFunctions.ese_score_index,
      prefix = prefix,
      genome_build = genome_build,
      force_rerun_token = force_rerun_token,
      docker = introme_functions_docker,
      preemptible = preemptible
  }
  call Tasks.EnsembleInfer {
    input:
      splicing_anno_vcf = SplicingAnno.splicing_anno_vcf,
      assets_tar = FetchIntromeAssets.assets_tar,
      prefix = prefix,
      commit_sha = FetchIntromeAssets.commit_sha,
      force_rerun_token = force_rerun_token,
      docker = introme_functions_docker,
      preemptible = heavy_preemptible
  }

  output {
    File introme_tsv = EnsembleInfer.introme_tsv
    File introme_vcf = EnsembleInfer.introme_vcf
    File introme_vcf_index = EnsembleInfer.introme_vcf_index
    String commit_sha = FetchIntromeAssets.commit_sha
  }
}
