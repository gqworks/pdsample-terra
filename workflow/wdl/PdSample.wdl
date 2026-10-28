version 1.0

import "tasks/SampleTasks.wdl" as Tasks
import "Introme.wdl" as IntromeWf
import "https://raw.githubusercontent.com/broadinstitute/gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl" as Mitochondria
import "vendor/gatk-sv/GATKSVPipelineSingleSample.wdl" as GatkSv

workflow PdSample {
  meta {
    description: "One Terra sample row. Every call is made from a BAM (+ .bai). GATK-SV v1.1.1 (Scramble on, MELT off) runs unless workspace sv_joint_vcf contains sample_id, in which case that joint-called SV VCF is subset and GATKSVPipelineSingleSample is skipped; if sample_id is absent from the joint VCF, GATK-SV still runs. A cohort is a Terra sample_set that launches this workflow once per row. CRAM inputs (with or without an index) are converted to BAM/BAI by EnsureBam."
  }
  input {
    String sample_id
    # Terra sample column bam_or_cram_file. May be BAM or CRAM (magic-byte detected).
    File bam
    # Optional. Required in practice for BAM when you already have a .bai; CRAM
    # without a .crai is fine — EnsureBam converts and indexes.
    File? bai
    String sex = ""

    File ref_fasta
    File ref_fasta_index
    File? ref_dict
    String requester_pays_project = ""

    File vep_cache
    String vep_extra_args = ""
    File? avi_table
    File? avi_table_index
    # Optional genome-wide CADD / REVEL plugins (fallback when not using dbNSFP 5.4a).
    File? cadd_snv
    File? cadd_snv_index
    File? cadd_indels
    File? cadd_indels_index
    File? revel
    File? revel_index
    # Academic dbNSFP 5.4a VEP-ready single file (bgzip + tabix), e.g. dbNSFP5.4a_grch38.gz.
    # Includes clinvar_clnsig + gnomAD4.1 joint AF (no separate ClinVar/gnomAD --custom).
    File? dbnsfp
    File? dbnsfp_index
    String dbnsfp_fields = "SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Eigen-phred_coding,MutationTaster_score,MutationTaster_pred,fathmm-XF_coding_score,fathmm-XF_coding_pred,CADD_phred,REVEL_score,clinvar_clnsig,gnomAD4.1_joint_AF,gnomAD4.1_joint_POPMAX_AF"
    # Off by default: SIFT/PolyPhen2 come from dbNSFP 5.4a (avoid duplicating cache-native VEP).
    Boolean vep_sift_polyphen = false

    File somalier_sites
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
    String mmsplice_cpu_docker = "headoncollusion/mmsplice:v2.1-cpu"
    String mmsplice_gpu_docker = "headoncollusion/mmsplice:v2.1-gpu"
    String pangolin_cpu_docker = "headoncollusion/pangolin:v1.3-cpu"
    String pangolin_gpu_docker = "headoncollusion/pangolin:v1.3-gpu"
    String spip_docker = "headoncollusion/spip:v1.1"
    String spliceogen_docker = "headoncollusion/spliceogen:v2.0"
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
    String? ethnicity
    File? phenotype_report

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
    String sv_scramble_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/scramble:2024-10-25-v0.29-beta-5ea22a52"
    String sv_wham_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/wham:2024-10-25-v0.29-beta-5ea22a52"
    String sv_sv_base_mini_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52"
    String sv_sv_base_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base:2024-10-25-v0.29-beta-5ea22a52"
    String sv_sv_pipeline_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-10-02-v1.1-483973d6"
    String sv_sv_pipeline_qc_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-10-02-v1.1-483973d6"
    String sv_linux_docker = "marketplace.gcr.io/google/ubuntu1804"
    String sv_cnmops_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/cnmops:2025-09-02-v1.0.5-f091af0b"
    String sv_gatk_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT"
    String sv_condense_counts_docker = "us.gcr.io/broad-dsde-methods/tsharpe/gatk:4.2.6.1-57-g9e03432"
    String sv_genomes_in_the_cloud_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"
    String sv_samtools_cloud_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/samtools-cloud:2024-10-25-v0.29-beta-5ea22a52"
    String sv_cloud_sdk_docker = "google/cloud-sdk"
    String sv_gq_recalibrator_gatk_docker = "us.gcr.io/broad-dsde-methods/markw/gatk:mw-tb-form-sv-filter-training-data-eb62a3"
    # v1.1.1 defaults: Scramble on, MELT off (MELT is deprecated).
    Boolean sv_use_melt = false
    Boolean sv_use_scramble = true
    # Optional GatherSampleEvidence outputs for this sample. When set, GATK-SV
    # skips the matching caller / coverage / PE-SR collection. Provide all of
    # counts + pe + sr + sd to skip evidence collection entirely; provide each
    # caller VCF to skip that caller. Indexes are not separate WDL inputs —
    # keep .tbi objects next to the primary paths in GCS when tools need them.
    File? sv_case_manta_vcf
    File? sv_case_wham_vcf
    File? sv_case_scramble_vcf
    File? sv_case_counts_file
    File? sv_case_pe_file
    File? sv_case_sr_file
    File? sv_case_sd_file
    # Optional cohort joint-called SV VCF (e.g. GATK-SV CleanVcf / FinalizeGVCF
    # multi-sample output). Bind as a workspace attribute (workspace.sv_joint_vcf),
    # not a sample-table column. When set and sample_id is in the header,
    # SubsetJointSvVcf keeps only that sample and GATKSVPipelineSingleSample is
    # skipped. If sample_id is absent from the joint VCF, fall back to
    # GATKSVPipelineSingleSample. Panel / evidence inputs below are unused on
    # the joint path but remain required File inputs for the single-sample path.
    File? sv_joint_vcf
    File? sv_joint_vcf_index
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
    # Manta runs one local job per core and accounts 2 GiB to each job.
    # 16 cores with jobs_per_cpu 1.0 is 16 jobs, so the VM is 32 GiB (2 GiB
    # per vCPU). Upstream default is 8 cores. Wham -x and samtools -@ read
    # the same core count; Wham's default is 4.
    Float sv_manta_jobs_per_cpu = 1.0
    Int sv_caller_cpu_cores = 16
    Float sv_caller_mem_gb = 32.0
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
    File sv_ref_std_scramble_vcf_tar
    File sv_ref_std_wham_vcf_tar
    File sv_ref_panel_del_bed
    File sv_ref_panel_dup_bed
    File sv_depth_exclude_list
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
    File sv_clustering_config_part1
    File sv_clustering_config_part2
    File sv_stratification_config_part1
    File sv_stratification_config_part2
    File sv_gq_recalibrator_model_file
    Array[File] sv_clustering_track_bed_files = [
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.SimpRep.sorted.pad_100.merged.bed",
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.SegDup.sorted.merged.bed",
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.RM.sorted.merged.bed"
    ]
    Array[String] sv_clustering_track_names = ["SR", "SD", "RM"]
    Array[File] sv_genome_tracks = [
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-RepeatMasker.bed.gz",
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Segmental-Dups.bed.gz",
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38-Simple-Repeats.bed.gz",
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s100.bed.gz",
      "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/ucsc-genome-tracks/hg38_umap_s24.bed.gz"
    ]
    Array[String] sv_recalibrate_gq_args = [
      "--keep-homvar false",
      "--keep-homref true",
      "--keep-multiallelic true",
      "--skip-genotype-filtering true",
      "--min-samples-to-estimate-allele-frequency -1"
    ]
    String sv_sl_filter_args = "--gq-scale-factor 1.012072 --upper-sl-cap 1500 --lower-sl-cap -1000 --sl-shift 300"
    File sv_hervk_reference = "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/HERVK.sorted.bed.gz"
    File sv_line1_reference = "gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/LINE1.sorted.bed.gz"
    Int sv_refine_complex_n_per_split = 15000

    Int preemptible = 3
    String hc_gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    String python_docker = "python:3.11.11-slim-bookworm"
  }

  call Tasks.FetchRepoAssets {
    input:
      git_commit = pdvar_git_commit,
      python_docker = python_docker
  }

  scatter (unused in if defined(ref_dict) then [1] else []) {
    call Tasks.CopyReferenceDict as ProvidedReferenceDict {
      input:
        ref_dict = select_first([ref_dict]),
        preemptible = preemptible
    }
  }
  scatter (unused in if !defined(ref_dict) then [1] else []) {
    call Tasks.SamtoolsReferenceDict as BuiltReferenceDict {
      input:
        ref_fasta = ref_fasta,
        preemptible = preemptible
    }
  }
  Array[File] reference_dict_files = flatten([ProvidedReferenceDict.dict, BuiltReferenceDict.dict])
  # Prefer the workflow input path so Cromwell co-localizes Homo_sapiens_assembly38.dict
  # beside the FASTA. CopyReferenceDict renames to reference.dict, which breaks imported
  # GATK tasks (mito, and HC with plain -R) that look for <fasta_basename>.dict next door.
  # GenotypeHardFilter still uses reference_dict_files and bundles in-task.
  File gatk_ref_dict = select_first([ref_dict, select_first(reference_dict_files)])

  # Normalize BAM or CRAM (+ optional index) to sample_id.bam / sample_id.bam.bai
  # before any reader. CRAM without a .crai is converted and indexed here.
  call Tasks.EnsureBam {
    input:
      alignment = bam,
      alignment_index = bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sample_id = sample_id,
      preemptible = preemptible
  }
  File sample_bam = EnsureBam.bam
  File sample_bai = EnsureBam.bai

  call Tasks.ExtractBamHeader {
    input:
      bam = sample_bam,
      bam_index = sample_bai,
      sample_id = sample_id,
      preemptible = preemptible
  }
  call Tasks.CheckReference {
    input:
      bam_header = ExtractBamHeader.header,
      reference_dict = select_first(reference_dict_files),
      check_reference_py = FetchRepoAssets.check_reference_py,
      sample_id = sample_id,
      python_docker = python_docker
  }
  call Tasks.VerifyBamID {
    input:
      bam = sample_bam,
      bam_index = sample_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      svd_ud = FetchRepoAssets.verifybamid_ud,
      svd_v = FetchRepoAssets.verifybamid_v,
      svd_bed = FetchRepoAssets.verifybamid_bed,
      svd_mu = FetchRepoAssets.verifybamid_mu,
      sample_id = sample_id,
      reference_report = CheckReference.report,
      preemptible = preemptible
  }
  call Tasks.SomalierCheck {
    input:
      bam = sample_bam,
      bam_index = sample_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sites = somalier_sites,
      sample_id = sample_id,
      reference_report = CheckReference.report,
      preemptible = preemptible
  }
  call Tasks.Mosdepth {
    input:
      bam = sample_bam,
      bam_index = sample_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      bed = FetchRepoAssets.panel_bed,
      sample_id = sample_id,
      reference_report = CheckReference.report,
      preemptible = preemptible
  }
  call Tasks.QcGate {
    input:
      qc_summary_py = FetchRepoAssets.qc_summary_py,
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
      bam = sample_bam,
      bam_index = sample_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = gatk_ref_dict,
      interval_bed = FetchRepoAssets.panel_bed,
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
      ref_dict = select_first(reference_dict_files),
      interval_bed = FetchRepoAssets.panel_bed,
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
  call Tasks.AnnotateAvi {
    input:
      vcf = NormPass.normalized_vcf,
      annotate_avi_py = FetchRepoAssets.annotate_avi_py,
      avi_table = avi_table,
      avi_table_index = avi_table_index,
      preemptible = preemptible
  }
  call Tasks.VepAnnotate {
    input:
      vcf = NormPass.normalized_vcf,
      vcf_index = NormPass.normalized_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      cache_dir = vep_cache,
      cadd_snv = cadd_snv,
      cadd_snv_index = cadd_snv_index,
      cadd_indels = cadd_indels,
      cadd_indels_index = cadd_indels_index,
      revel = revel,
      revel_index = revel_index,
      dbnsfp = dbnsfp,
      dbnsfp_index = dbnsfp_index,
      dbnsfp_fields = dbnsfp_fields,
      vep_sift_polyphen = vep_sift_polyphen,
      extra_args = vep_extra_args,
      preemptible = preemptible
  }
  call Tasks.LiftGnomadInfo {
    input:
      vcf = VepAnnotate.vep_vcf,
      vcf_index = VepAnnotate.vep_index,
      lift_gnomad_from_csq_py = FetchRepoAssets.lift_gnomad_from_csq_py,
      sample_id = sample_id,
      preemptible = preemptible
  }
  call Tasks.ExpansionHunter {
    input:
      bam = sample_bam,
      bam_index = sample_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      catalog = FetchRepoAssets.eh_catalog,
      sample_id = sample_id,
      sex = QcGate.sex,
      qc_gate = QcGate.gate,
      preemptible = preemptible
  }
  call Tasks.Gauchian {
    input:
      bam = sample_bam,
      bam_index = sample_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      sample_id = sample_id,
      qc_gate = QcGate.gate,
      gauchian_docker = gauchian_docker,
      preemptible = preemptible
  }
  # Gate mito / GATK-SV on QC. EnsureBam already emits sample_id.bam / .bam.bai;
  # ReleaseAfterQc hardlinks again so Manta 1.6 does not start before the gate.
  call Tasks.ReleaseAfterQc as ReleaseBam {
    input:
      payload = sample_bam,
      gate = QcGate.gate,
      name = sample_id + ".bam"
  }
  call Tasks.ReleaseAfterQc as ReleaseBai {
    input:
      payload = sample_bai,
      gate = QcGate.gate,
      name = sample_id + ".bam.bai"
  }
  call Mitochondria.MitochondriaPipeline as Mito {
    input:
      wgs_aligned_input_bam_or_cram = ReleaseBam.released,
      wgs_aligned_input_bam_or_cram_index = ReleaseBai.released,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = gatk_ref_dict,
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
      # SplitMultiAllelicsAndRemoveNonPassSites reads InitialFilter's VCF and
      # does not take its index. GATK rejects that file when it is block-compressed
      # ("drivingVariantFile ... Try running IndexFeatureFile"). Leave the
      # upstream default (uncompressed .vcf plus .idx).
      compress_output_vcf = false,
      preemptible_tries = preemptible
  }
  # When sv_joint_vcf is set, probe/subset first. If sample_id is in the header,
  # use that subset and skip GatkSvSingle. If absent, fall back to GatkSvSingle.
  # scatter+if keeps Array[File] outputs (same pattern as reference_dict_files).
  scatter (unused in if defined(sv_joint_vcf) then [1] else []) {
    call Tasks.SubsetJointSvVcf {
      input:
        joint_vcf = select_first([sv_joint_vcf]),
        joint_vcf_index = sv_joint_vcf_index,
        sample_id = sample_id,
        preemptible = preemptible
    }
  }
  Boolean use_joint_sv = select_first(flatten([SubsetJointSvVcf.sample_found, [false]]))
  scatter (unused in if !use_joint_sv then [1] else []) {
    call GatkSv.GATKSVPipelineSingleSample as GatkSvSingle {
      input:
        sample_id = sample_id,
        batch = "pd-sample",
        bam_or_cram_file = ReleaseBam.released,
        bam_or_cram_index = ReleaseBai.released,
        reference_version = "38",
        use_melt = sv_use_melt,
        use_scramble = sv_use_scramble,
        case_manta_vcf = sv_case_manta_vcf,
        case_wham_vcf = sv_case_wham_vcf,
        case_scramble_vcf = sv_case_scramble_vcf,
        case_counts_file = sv_case_counts_file,
        case_pe_file = sv_case_pe_file,
        case_sr_file = sv_case_sr_file,
        case_sd_file = sv_case_sd_file,
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
        gq_recalibrator_gatk_docker = sv_gq_recalibrator_gatk_docker,
        preprocessed_intervals = sv_preprocessed_intervals,
        manta_region_bed = sv_manta_region_bed,
        manta_region_bed_index = sv_manta_region_bed_index,
        manta_jobs_per_cpu = sv_manta_jobs_per_cpu,
        runtime_attr_manta = object {
          cpu_cores: sv_caller_cpu_cores,
          mem_gb: sv_caller_mem_gb
        },
        runtime_attr_wham = object {
          cpu_cores: sv_caller_cpu_cores,
          mem_gb: sv_caller_mem_gb
        },
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
        ref_std_scramble_vcf_tar = sv_ref_std_scramble_vcf_tar,
        ref_std_wham_vcf_tar = sv_ref_std_wham_vcf_tar,
        ref_panel_del_bed = sv_ref_panel_del_bed,
        ref_panel_dup_bed = sv_ref_panel_dup_bed,
        depth_exclude_list = sv_depth_exclude_list,
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
        clustering_config_part1 = sv_clustering_config_part1,
        clustering_config_part2 = sv_clustering_config_part2,
        stratification_config_part1 = sv_stratification_config_part1,
        stratification_config_part2 = sv_stratification_config_part2,
        clustering_track_bed_files = sv_clustering_track_bed_files,
        clustering_track_names = sv_clustering_track_names,
        gq_recalibrator_model_file = sv_gq_recalibrator_model_file,
        genome_tracks = sv_genome_tracks,
        recalibrate_gq_args = sv_recalibrate_gq_args,
        sl_filter_args = sv_sl_filter_args,
        HERVK_reference = sv_hervk_reference,
        LINE1_reference = sv_line1_reference,
        refine_complex_n_per_split = sv_refine_complex_n_per_split,
        svtk_to_gatk_script = FetchRepoAssets.svtk_to_gatk_script,
        manta_docker = sv_manta_docker,
        scramble_docker = sv_scramble_docker,
        wham_docker = sv_wham_docker
    }
  }
  # Prefer GatkSv when both ran (sample absent from joint → fallback). Prefer
  # the joint subset when GatkSv was skipped (sample_found).
  File sample_sv_vcf = select_first(flatten([GatkSvSingle.final_vcf, SubsetJointSvVcf.vcf]))
  File sample_sv_vcf_index = select_first(flatten([GatkSvSingle.final_vcf_idx, SubsetJointSvVcf.vcf_index]))
  # Stage as sample_id.sv.vcf.gz (GATK-SV finals are opaque; joint subset already matches).
  call Tasks.ReleaseAfterQc as NameSvVcf {
    input:
      payload = sample_sv_vcf,
      gate = QcGate.gate,
      name = sample_id + ".sv.vcf.gz"
  }
  call Tasks.ReleaseAfterQc as NameSvVcfIndex {
    input:
      payload = sample_sv_vcf_index,
      gate = QcGate.gate,
      name = sample_id + ".sv.vcf.gz.tbi"
  }
  call IntromeWf.Introme as Introme {
    input:
      vcf = LiftGnomadInfo.lifted_vcf,
      vcf_index = LiftGnomadInfo.lifted_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gtf = introme_gtf,
      bed = FetchRepoAssets.panel_bed,
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
      sv_vcf = NameSvVcf.released,
      sv_status = if use_joint_sv then "joint" else "called",
      expansionhunter_vcf = ExpansionHunter.vcf,
      gauchian_json = Gauchian.json,
      mtdna_vcf = Mito.out_vcf,
      introme_tsv = Introme.introme_tsv,
      avi_tsv = AnnotateAvi.avi_tsv,
      qc_json = QcGate.qc_json,
      config_yaml = FetchRepoAssets.config_yaml,
      columns_yaml = FetchRepoAssets.columns_yaml,
      phenotypes_yaml = FetchRepoAssets.phenotypes_yaml,
      panel_csv = FetchRepoAssets.panel_csv,
      gene_bed = FetchRepoAssets.gene_bed,
      exon_bed = FetchRepoAssets.exon_bed,
      expression_tsv = expression_tsv,
      splicevault_tsv = splicevault_tsv,
      prs_weights = prs_weights,
      age_onset = age_onset,
      symptoms = symptoms,
      ethnicity = ethnicity,
      phenotype_report = phenotype_report,
      pdvar_docker = pdvar_docker,
      pdvar_git_url = pdvar_git_url,
      pdvar_git_commit = pdvar_git_commit
  }

  output {
    Array[File] candidates = Assess.candidates
    File vep_vcf = LiftGnomadInfo.lifted_vcf
    File vep_vcf_index = LiftGnomadInfo.lifted_vcf_index
    File sv_vcf = NameSvVcf.released
    File sv_vcf_index = NameSvVcfIndex.released
  }
}
