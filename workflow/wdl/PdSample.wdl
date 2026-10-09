version 1.0

import "tasks/SampleTasks.wdl" as Tasks
import "Introme.wdl" as IntromeWf
import "https://raw.githubusercontent.com/broadinstitute/gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl" as Mitochondria
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.29-beta/wdl/GATKSVPipelineSingleSample.wdl" as GatkSv

workflow PdSample {
  meta {
    description: "One Terra sample row. Every call is made from the CRAM, including required GATK-SV v0.29-beta. A cohort is a Terra sample_set that launches this workflow once per row."
  }
  input {
    String sample_id
    File cram
    File crai
    String sex = ""

    File ref_fasta
    File ref_fasta_index
    File? ref_dict
    String requester_pays_project = ""

    File panel_bed
    File panel_csv
    File gene_bed
    File exon_bed
    File eh_catalog
    File config_yaml
    File columns_yaml
    File phenotypes_yaml
    File check_reference_py
    File qc_summary_py

    File vep_cache
    String vep_extra_args = ""
    File? gnomad_vcf
    File? gnomad_vcf_index
    File? clinvar_vcf
    File? clinvar_vcf_index

    File somalier_sites
    File verifybamid_ud
    File verifybamid_v
    File verifybamid_bed
    File verifybamid_mu
    Float max_freemix = 0.02
    Float min_mean_coverage = 0
    Boolean hard_fail_contamination = true
    Boolean hard_fail_sex = false
    Boolean hard_fail_coverage = false

    Boolean wgs_small_variants = false

    File introme_gtf
    File? introme_assets_tar
    Boolean introme_quality_filter = false
    Int introme_min_qual = 200
    Int introme_min_dp = 20
    Int introme_min_ad = 5
    Boolean introme_filter_rare = true
    Float introme_allele_frequency = 0.01
    String introme_af_field = ""
    Boolean introme_restrict_to_panel = true
    Boolean introme_use_gpu = false
    String introme_gpu_type = "nvidia-tesla-t4"
    Int introme_gpu_count = 1
    String introme_nvidia_driver_version = "535.183.01"
    File? introme_spliceai_db
    File? introme_pangolin_db
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

    String gauchian_docker = "python:3.11.11-slim-bookworm"
    String pdvar_docker = "python:3.11.11-slim-bookworm"
    String pdvar_git_url = "https://github.com/gqworks/pdsample-terra.git"
    String pdvar_git_commit = "main"
    File? expression_tsv
    File? splicevault_tsv
    File? prs_weights
    String? age_onset
    String? symptoms

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa
    File blacklisted_sites
    File blacklisted_sites_index
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa
    File shift_back_chain
    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list

    String sv_manta_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/manta:2023-09-14-v0.28.3-beta-3f22f94d"
    String sv_scramble_docker = "us.gcr.io/broad-dsde-methods/markw/scramble:mw-scramble-99af4c50"
    String sv_wham_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/wham:2024-01-24-v0.28.4-beta-9debd6d7"
    String sv_sv_base_mini_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-01-24-v0.28.4-beta-9debd6d7"
    String sv_sv_base_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base:2024-01-24-v0.28.4-beta-9debd6d7"
    String sv_sv_pipeline_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-08-19-v0.28.5-beta-84a0627d"
    String sv_sv_pipeline_qc_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-08-19-v0.28.5-beta-84a0627d"
    String sv_linux_docker = "marketplace.gcr.io/google/ubuntu1804"
    String sv_cnmops_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/cnmops:2024-06-04-v0.28.5-beta-a8dfecba"
    String sv_gatk_docker = "us.gcr.io/broad-dsde-methods/eph/gatk:2024-07-02-4.6.0.0-1-g4af2b49e9-NIGHTLY-SNAPSHOT"
    String sv_condense_counts_docker = "us.gcr.io/broad-dsde-methods/tsharpe/gatk:4.2.6.1-57-g9e03432"
    String sv_genomes_in_the_cloud_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
    String sv_samtools_cloud_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/samtools-cloud:2024-01-24-v0.28.4-beta-9debd6d7"
    String sv_cloud_sdk_docker = "google/cloud-sdk"
    File sv_ref_samples_list
    File sv_ref_ped_file
    File sv_genome_file
    File sv_primary_contigs_list
    File sv_primary_contigs_fai
    File sv_reference_fasta
    File sv_reference_index
    File sv_reference_dict
    File sv_ref_panel_vcf
    File sv_autosome_file
    File sv_allosome_file
    File sv_preprocessed_intervals
    File sv_manta_region_bed
    File sv_manta_region_bed_index
    File sv_sd_locs_vcf
    File sv_wham_include_list_bed_file
    Boolean sv_run_vcf_qc = false
    File sv_wgd_scoring_mask
    Int sv_min_svsize = 50
    File sv_contig_ploidy_model_tar
    File sv_gcnv_model_tars_list
    File sv_ref_panel_bincov_matrix
    File sv_ref_pesr_disc_files_list
    File sv_ref_pesr_split_files_list
    File sv_ref_pesr_sd_files_list
    Int sv_ref_copy_number_autosomal_contigs = 2
    Int sv_gcnv_qs_cutoff = 30
    File sv_cnmops_exclude_list
    Int sv_matrix_qc_distance = 1000000
    File sv_ref_std_manta_vcf_tar
    File sv_ref_std_wham_vcf_tar
    File sv_ref_panel_del_bed
    File sv_ref_panel_dup_bed
    File sv_depth_exclude_list
    File sv_empty_file
    Float sv_depth_exclude_overlap_fraction = 0.5
    Float sv_depth_interval_overlap = 0.8
    File sv_pesr_exclude_intervals
    Float sv_pesr_interval_overlap = 0.1
    Int sv_pesr_breakend_window = 300
    File sv_rmsk
    File sv_segdups
    Int sv_genotyping_n_per_split = 1000
    Int sv_n_RD_genotype_bins = 100000
    File sv_cutoffs
    File sv_genotype_pesr_pesr_sepcutoff
    File sv_genotype_pesr_depth_sepcutoff
    File sv_genotype_depth_pesr_sepcutoff
    File sv_genotype_depth_depth_sepcutoff
    File sv_SR_metrics
    File sv_PE_metrics
    File sv_bin_exclude
    Float sv_clean_vcf_min_sr_background_fail_batches = 0.5
    File sv_cytobands
    File sv_mei_bed
    Int sv_max_shard_size_resolve = 500
    Int sv_clean_vcf_max_shards_per_chrom_clean_vcf_step1 = 200
    Int sv_clean_vcf_min_records_per_shard_clean_vcf_step1 = 5000
    Int sv_clean_vcf_samples_per_clean_vcf_step2_shard = 100
    Int sv_clean_vcf5_records_per_shard = 5000
    Int sv_clean_vcf1b_records_per_shard = 10000
    File sv_protein_coding_gtf
    File sv_noncoding_bed
    Int sv_annotation_sv_per_shard = 5000
    File sv_qc_definitions

    Int preemptible = 3
    String hc_gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    String python_docker = "python:3.11.11-slim-bookworm"
  }

  call Tasks.MakeReferenceDict {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      preemptible = preemptible
  }
  call Tasks.ExtractCramHeader {
    input:
      cram = cram,
      cram_index = crai,
      sample_id = sample_id,
      preemptible = preemptible
  }
  call Tasks.CheckReference {
    input:
      cram_header = ExtractCramHeader.header,
      reference_dict = MakeReferenceDict.dict,
      check_reference_py = check_reference_py,
      sample_id = sample_id,
      python_docker = python_docker
  }
  call Tasks.VerifyBamID {
    input:
      cram = cram,
      cram_index = crai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      svd_ud = verifybamid_ud,
      svd_v = verifybamid_v,
      svd_bed = verifybamid_bed,
      svd_mu = verifybamid_mu,
      sample_id = sample_id,
      reference_report = CheckReference.report,
      preemptible = preemptible
  }
  call Tasks.SomalierCheck {
    input:
      cram = cram,
      cram_index = crai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sites = somalier_sites,
      sample_id = sample_id,
      reference_report = CheckReference.report,
      preemptible = preemptible
  }
  call Tasks.Mosdepth {
    input:
      cram = cram,
      cram_index = crai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      bed = panel_bed,
      sample_id = sample_id,
      reference_report = CheckReference.report,
      preemptible = preemptible
  }
  call Tasks.QcGate {
    input:
      qc_summary_py = qc_summary_py,
      selfsm = VerifyBamID.self_sm,
      somalier_samples = SomalierCheck.samples_tsv,
      mosdepth_summary = Mosdepth.summary,
      reference_report = CheckReference.report,
      sample_id = sample_id,
      provided_sex = sex,
      max_freemix = max_freemix,
      min_mean_coverage = min_mean_coverage,
      hard_fail_contamination = hard_fail_contamination,
      hard_fail_sex = hard_fail_sex,
      hard_fail_coverage = hard_fail_coverage,
      python_docker = python_docker
  }
  call Tasks.HaplotypeCaller {
    input:
      cram = cram,
      cram_index = crai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = MakeReferenceDict.dict,
      interval_bed = panel_bed,
      use_intervals = !wgs_small_variants,
      sample_id = sample_id,
      qc_gate = QcGate.gate,
      requester_pays_project = requester_pays_project,
      gatk_docker = hc_gatk_docker,
      preemptible = preemptible
  }
  call Tasks.GenotypeHardFilter {
    input:
      gvcf = HaplotypeCaller.gvcf,
      gvcf_index = HaplotypeCaller.gvcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = MakeReferenceDict.dict,
      interval_bed = panel_bed,
      use_intervals = !wgs_small_variants,
      sample_id = sample_id,
      gatk_docker = hc_gatk_docker,
      preemptible = preemptible
  }
  call Tasks.NormPass {
    input:
      vcf = GenotypeHardFilter.vcf,
      vcf_index = GenotypeHardFilter.vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sample_id = sample_id,
      preemptible = preemptible
  }
  call Tasks.VepAnnotate {
    input:
      vcf = NormPass.normalized_vcf,
      vcf_index = NormPass.normalized_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      cache_dir = vep_cache,
      gnomad_vcf = gnomad_vcf,
      gnomad_vcf_index = gnomad_vcf_index,
      clinvar_vcf = clinvar_vcf,
      clinvar_vcf_index = clinvar_vcf_index,
      extra_args = vep_extra_args,
      preemptible = preemptible
  }
  call Tasks.LiftGnomadInfo {
    input:
      vcf = VepAnnotate.vep_vcf,
      vcf_index = VepAnnotate.vep_index,
      preemptible = preemptible
  }
  call Tasks.ExpansionHunter {
    input:
      cram = cram,
      cram_index = crai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      catalog = eh_catalog,
      sample_id = sample_id,
      sex = QcGate.sex,
      qc_gate = QcGate.gate,
      preemptible = preemptible
  }
  call Tasks.Gauchian {
    input:
      cram = cram,
      cram_index = crai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sample_id = sample_id,
      qc_gate = QcGate.gate,
      gauchian_docker = gauchian_docker,
      preemptible = preemptible
  }
  call Tasks.ReleaseAfterQc as ReleaseCram {
    input:
      payload = cram,
      gate = QcGate.gate,
      name = sample_id + ".cram"
  }
  call Tasks.ReleaseAfterQc as ReleaseCrai {
    input:
      payload = crai,
      gate = QcGate.gate,
      name = sample_id + ".cram.crai"
  }
  call Mitochondria.MitochondriaPipeline as Mito {
    input:
      wgs_aligned_input_bam_or_cram = ReleaseCram.released,
      wgs_aligned_input_bam_or_cram_index = ReleaseCrai.released,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = MakeReferenceDict.dict,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      mt_shifted_dict = mt_shifted_dict,
      mt_shifted_fasta = mt_shifted_fasta,
      mt_shifted_fasta_index = mt_shifted_fasta_index,
      mt_shifted_amb = mt_shifted_amb,
      mt_shifted_ann = mt_shifted_ann,
      mt_shifted_bwt = mt_shifted_bwt,
      mt_shifted_pac = mt_shifted_pac,
      mt_shifted_sa = mt_shifted_sa,
      shift_back_chain = shift_back_chain,
      control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
      non_control_region_interval_list = non_control_region_interval_list,
      requester_pays_project = requester_pays_project,
      compress_output_vcf = true,
      preemptible_tries = preemptible
  }
  call GatkSv.GATKSVPipelineSingleSample as GatkSvSingle {
    input:
      sample_id = sample_id,
      batch = "pd-sample",
      bam_or_cram_file = ReleaseCram.released,
      bam_or_cram_index = ReleaseCrai.released,
      reference_version = "38",
      ref_samples_list = sv_ref_samples_list,
      ref_ped_file = sv_ref_ped_file,
      genome_file = sv_genome_file,
      primary_contigs_list = sv_primary_contigs_list,
      primary_contigs_fai = sv_primary_contigs_fai,
      reference_fasta = sv_reference_fasta,
      reference_index = sv_reference_index,
      reference_dict = sv_reference_dict,
      ref_panel_vcf = sv_ref_panel_vcf,
      autosome_file = sv_autosome_file,
      allosome_file = sv_allosome_file,
      sv_base_mini_docker = sv_sv_base_mini_docker,
      sv_base_docker = sv_sv_base_docker,
      sv_pipeline_docker = sv_sv_pipeline_docker,
      sv_pipeline_qc_docker = sv_sv_pipeline_qc_docker,
      linux_docker = sv_linux_docker,
      cnmops_docker = sv_cnmops_docker,
      gatk_docker = sv_gatk_docker,
      condense_counts_docker = sv_condense_counts_docker,
      genomes_in_the_cloud_docker = sv_genomes_in_the_cloud_docker,
      samtools_cloud_docker = sv_samtools_cloud_docker,
      cloud_sdk_docker = sv_cloud_sdk_docker,
      preprocessed_intervals = sv_preprocessed_intervals,
      manta_region_bed = sv_manta_region_bed,
      manta_region_bed_index = sv_manta_region_bed_index,
      sd_locs_vcf = sv_sd_locs_vcf,
      wham_include_list_bed_file = sv_wham_include_list_bed_file,
      run_vcf_qc = sv_run_vcf_qc,
      wgd_scoring_mask = sv_wgd_scoring_mask,
      min_svsize = sv_min_svsize,
      contig_ploidy_model_tar = sv_contig_ploidy_model_tar,
      gcnv_model_tars_list = sv_gcnv_model_tars_list,
      ref_panel_bincov_matrix = sv_ref_panel_bincov_matrix,
      ref_pesr_disc_files_list = sv_ref_pesr_disc_files_list,
      ref_pesr_split_files_list = sv_ref_pesr_split_files_list,
      ref_pesr_sd_files_list = sv_ref_pesr_sd_files_list,
      ref_copy_number_autosomal_contigs = sv_ref_copy_number_autosomal_contigs,
      gcnv_qs_cutoff = sv_gcnv_qs_cutoff,
      cnmops_exclude_list = sv_cnmops_exclude_list,
      matrix_qc_distance = sv_matrix_qc_distance,
      ref_std_manta_vcf_tar = sv_ref_std_manta_vcf_tar,
      ref_std_wham_vcf_tar = sv_ref_std_wham_vcf_tar,
      ref_panel_del_bed = sv_ref_panel_del_bed,
      ref_panel_dup_bed = sv_ref_panel_dup_bed,
      depth_exclude_list = sv_depth_exclude_list,
      empty_file = sv_empty_file,
      depth_exclude_overlap_fraction = sv_depth_exclude_overlap_fraction,
      depth_interval_overlap = sv_depth_interval_overlap,
      pesr_exclude_intervals = sv_pesr_exclude_intervals,
      pesr_interval_overlap = sv_pesr_interval_overlap,
      pesr_breakend_window = sv_pesr_breakend_window,
      rmsk = sv_rmsk,
      segdups = sv_segdups,
      genotyping_n_per_split = sv_genotyping_n_per_split,
      n_RD_genotype_bins = sv_n_RD_genotype_bins,
      cutoffs = sv_cutoffs,
      genotype_pesr_pesr_sepcutoff = sv_genotype_pesr_pesr_sepcutoff,
      genotype_pesr_depth_sepcutoff = sv_genotype_pesr_depth_sepcutoff,
      genotype_depth_pesr_sepcutoff = sv_genotype_depth_pesr_sepcutoff,
      genotype_depth_depth_sepcutoff = sv_genotype_depth_depth_sepcutoff,
      SR_metrics = sv_SR_metrics,
      PE_metrics = sv_PE_metrics,
      bin_exclude = sv_bin_exclude,
      clean_vcf_min_sr_background_fail_batches = sv_clean_vcf_min_sr_background_fail_batches,
      cytobands = sv_cytobands,
      mei_bed = sv_mei_bed,
      max_shard_size_resolve = sv_max_shard_size_resolve,
      clean_vcf_max_shards_per_chrom_clean_vcf_step1 = sv_clean_vcf_max_shards_per_chrom_clean_vcf_step1,
      clean_vcf_min_records_per_shard_clean_vcf_step1 = sv_clean_vcf_min_records_per_shard_clean_vcf_step1,
      clean_vcf_samples_per_clean_vcf_step2_shard = sv_clean_vcf_samples_per_clean_vcf_step2_shard,
      clean_vcf5_records_per_shard = sv_clean_vcf5_records_per_shard,
      clean_vcf1b_records_per_shard = sv_clean_vcf1b_records_per_shard,
      protein_coding_gtf = sv_protein_coding_gtf,
      noncoding_bed = sv_noncoding_bed,
      annotation_sv_per_shard = sv_annotation_sv_per_shard,
      qc_definitions = sv_qc_definitions,
      manta_docker = sv_manta_docker,
      scramble_docker = sv_scramble_docker,
      wham_docker = sv_wham_docker
    }
  call IntromeWf.Introme as Introme {
    input:
      vcf = LiftGnomadInfo.lifted_vcf,
      vcf_index = LiftGnomadInfo.lifted_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gtf = introme_gtf,
      bed = panel_bed,
      introme_assets_tar = introme_assets_tar,
      prefix = sample_id,
      restrict_to_panel = introme_restrict_to_panel,
      quality_filter = introme_quality_filter,
      min_qual = introme_min_qual,
      min_dp = introme_min_dp,
      min_ad = introme_min_ad,
      filter_rare = introme_filter_rare,
      allele_frequency = introme_allele_frequency,
      af_field = introme_af_field,
      use_gpu = introme_use_gpu,
      gpu_type = introme_gpu_type,
      gpu_count = introme_gpu_count,
      nvidia_driver_version = introme_nvidia_driver_version,
      spliceai_db_file = introme_spliceai_db,
      pangolin_db_file = introme_pangolin_db,
      data_preprocessing_docker = data_preprocessing_docker,
      variant_info_docker = variant_info_docker,
      spliceai_cpu_docker = spliceai_cpu_docker,
      spliceai_gpu_docker = spliceai_gpu_docker,
      mmsplice_cpu_docker = mmsplice_cpu_docker,
      mmsplice_gpu_docker = mmsplice_gpu_docker,
      pangolin_cpu_docker = pangolin_cpu_docker,
      pangolin_gpu_docker = pangolin_gpu_docker,
      spip_docker = spip_docker,
      spliceogen_docker = spliceogen_docker,
      introme_functions_docker = introme_functions_docker,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.Assess {
    input:
      sample_id = sample_id,
      sex = QcGate.sex,
      vep_vcf = LiftGnomadInfo.lifted_vcf,
      small_vcf = NormPass.normalized_vcf,
      sv_vcf = GatkSvSingle.final_vcf,
      expansionhunter_vcf = ExpansionHunter.vcf,
      gauchian_json = Gauchian.json,
      mtdna_vcf = Mito.out_vcf,
      introme_tsv = Introme.introme_tsv,
      qc_json = QcGate.qc_json,
      config_yaml = config_yaml,
      columns_yaml = columns_yaml,
      phenotypes_yaml = phenotypes_yaml,
      panel_csv = panel_csv,
      gene_bed = gene_bed,
      exon_bed = exon_bed,
      expression_tsv = expression_tsv,
      splicevault_tsv = splicevault_tsv,
      prs_weights = prs_weights,
      age_onset = age_onset,
      symptoms = symptoms,
      pdvar_docker = pdvar_docker,
      pdvar_git_url = pdvar_git_url,
      pdvar_git_commit = pdvar_git_commit
  }

  output {
    File candidates_tsv = Assess.candidates_tsv
    File report_html = Assess.report_html
    File report_md = Assess.report_md
    File qc_summary = Assess.qc_summary
    File liability_tsv = Assess.liability_tsv
    File followup_tsv = Assess.followup_tsv
    File gvcf = HaplotypeCaller.gvcf
    File gvcf_index = HaplotypeCaller.gvcf_index
    File small_vcf = NormPass.normalized_vcf
    File small_vcf_index = NormPass.normalized_vcf_index
    File vep_vcf = LiftGnomadInfo.lifted_vcf
    File vep_vcf_index = LiftGnomadInfo.lifted_vcf_index
    File sv_vcf = GatkSvSingle.final_vcf
    File sv_vcf_index = GatkSvSingle.final_vcf_idx
    String sv_status = "called"
    File expansionhunter_vcf = ExpansionHunter.vcf
    File expansionhunter_vcf_index = ExpansionHunter.vcf_index
    File gauchian_json = Gauchian.json
    File mtdna_vcf = Mito.out_vcf
    File mtdna_vcf_index = Mito.out_vcf_index
    File introme_tsv = Introme.introme_tsv
    File introme_vcf = Introme.introme_vcf
    String sex_used = QcGate.sex
    File reference_report = CheckReference.report
    File sample_json = Assess.sample_json
  }
}
