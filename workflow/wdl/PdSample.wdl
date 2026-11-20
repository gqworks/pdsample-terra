version 1.0

import "tasks/SampleTasks.wdl" as Tasks
import "Introme.wdl" as IntromeWf
import "vendor/mitochondria_m2_wdl/MitochondriaPipeline.wdl" as Mitochondria
import "vendor/gatk-sv/GATKSVPipelineSingleSample.wdl" as GatkSv

workflow PdSample {
  meta {
    description: "One Terra sample row. Every call is made from a BAM (+ .bai). Optional this.<name> priors (ProbePriorFile) skip only that file's producer: vep_vcf, small_vcf, sv_vcf, expansionhunter_vcf, gauchian_json, mtdna_vcf, introme, panel_csv, avi_tsv, gba1_pd_browser. Empty/placeholder/zero-byte priors run the caller. Else SV: sv_joint_vcf subset or GATK-SV v1.1.1 (Scramble on, MELT off). A cohort is a Terra sample_set that launches this workflow once per row. CRAM inputs (with or without an index) are converted to BAM/BAI by EnsureBam."
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
    # MitImpact 4 CSV for mtDNA ClinVar/gnomAD/MITOMAP/APOGEE2/MitoTIP annotation.
    # Unset leaves mtdna_annotation=not_assessed and passes the raw mito VCF through.
    File? mtdna_mitimpact
    String mtdna_mitimpact_release = "MitImpact4_V059"
    # Optional genome-wide CADD / REVEL plugins (fallback when not using dbNSFP 5.4a).
    File? cadd_snv
    File? cadd_snv_index
    File? cadd_indels
    File? cadd_indels_index
    File? revel
    File? revel_index
    # Academic dbNSFP 5.4a VEP-ready single file (bgzip + tabix), e.g. dbNSFP5.4a_grch38.gz.
    # Includes clinvar_clnsig + gnomAD4.1 joint AF for missense/splice; not genome-wide AF.
    File? dbnsfp
    File? dbnsfp_index
    String dbnsfp_fields = "SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Eigen-phred_coding,MutationTaster_score,MutationTaster_pred,fathmm-XF_coding_score,fathmm-XF_coding_pred,CADD_phred,REVEL_score,clinvar_clnsig,gnomAD4.1_joint_AF,gnomAD4.1_joint_POPMAX_AF"
    # Off by default: SIFT/PolyPhen2 come from dbNSFP 5.4a (avoid duplicating cache-native VEP).
    Boolean vep_sift_polyphen = false
    # Panel-sliced gnomAD v4.1 joint sites (AF_joint / AF_grpmax_joint). Rebuild when the
    # panel BED changes (workflow/scripts/build_panel_gnomad_af.sh). Unset keeps dbNSFP-only AF.
    File? gnomad_af_vcf
    File? gnomad_af_vcf_index
    # Panel-sliced ClinVar P/LP sites (CLNSIG Pathogenic / Likely_pathogenic).
    # Rebuild when the panel BED changes (workflow/scripts/build_panel_clinvar.sh).
    # Unset keeps dbNSFP missense/splice ClinVar only.
    File? clinvar_plp_vcf
    File? clinvar_plp_vcf_index

    File somalier_sites
    Float max_freemix = 0.02
    Float min_mean_coverage = 0
    Boolean hard_fail_contamination = true
    Boolean hard_fail_sex = false
    Boolean hard_fail_coverage = false

    Boolean wgs_small_variants = false

    File introme_gtf
    File? introme_assets_tar
    # Optional per-sample prior Introme TSV (Terra: this.introme). A real
    # non-empty Introme TSV skips the Introme subworkflow. Empty, missing,
    # placeholder, or zero-byte this.introme is treated as unset (Introme runs).
    # String so blank cells are not coerced to File before ProbePriorFile.
    String introme_prior = ""
    # Optional per-sample prior VEP VCF (Terra: this.vep_vcf / this.vep_vcf_index).
    # A real non-empty VCF skips VepAnnotate, AnnotateGnomadAf, LiftGnomadInfo,
    # and AnnotateClinvarSmall. Assess still runs. Empty/missing/placeholder/
    # zero-byte this.vep_vcf is treated as unset (VEP chain runs). Same probe
    # as introme_prior; WDL never size()s a bad File.
    String vep_prior_vcf = ""
    String vep_prior_vcf_index = ""
    # Same ProbePriorFile failsafe as introme_prior / vep_prior_vcf. Empty,
    # missing, placeholder, or zero-byte this.<name> runs the producer. A real
    # non-empty file skips only that producer. Assess still runs.
    String small_vcf_prior = ""
    String small_vcf_prior_index = ""
    String expansionhunter_vcf_prior = ""
    String gauchian_json_prior = ""
    String mtdna_vcf_prior = ""
    String panel_csv_prior = ""
    String avi_tsv_prior = ""
    String gba1_pd_browser_prior = ""
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
    # FetchRepoAssets is volatile and uses this ref for scripts/config/panel assets.
    # Default "main" so pushed script fixes are re-fetched; pin a SHA to freeze assets.
    String assets_git_commit = "main"
    # Pin a SHA so Assess call-cache invalidates when pdvar (src/pdvar) changes.
    # "main" does not change as an input string when the tip of main moves.
    String pdvar_git_commit = "d18d01efb2e7605d2e3e0f626e82f6b09fdf4219"
    File? prs_weights
    # Optional GP2-style phenotype text. Age at onset, ancestry, and clinical fields
    # are read from this file only (no separate sample-table phenotype columns).
    # String (not File?) so Terra blank sample cells ("") do not become an illegal
    # File path that fails Assess with "Could not build the path """.
    String phenotype_report = ""

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
    # Strings so Terra blank sample cells ("") stay unset instead of illegal File "".
    String sv_case_manta_vcf = ""
    String sv_case_wham_vcf = ""
    String sv_case_scramble_vcf = ""
    String sv_case_counts_file = ""
    String sv_case_pesr_disc = ""
    String sv_case_pesr_split = ""
    String sv_case_pesr_sd = ""
    # Optional per-sample prior SV VCF (Terra: this.sv_vcf / this.sv_vcf_index).
    # A real non-empty VCF skips SubsetJointSvVcf and GATKSVPipelineSingleSample
    # and reuses this callset as sv_vcf. Takes precedence over sv_joint_vcf.
    # Empty/missing/placeholder/zero-byte this.sv_vcf is probed like introme_prior.
    String sv_prior_vcf = ""
    String sv_prior_vcf_index = ""
    # Optional cohort joint-called SV VCF (e.g. GATK-SV CleanVcf / FinalizeGVCF
    # multi-sample output). Bind as a workspace attribute (workspace.sv_joint_vcf),
    # not a sample-table column. Used only when sv_prior_vcf is unset. When set
    # and sample_id is in the header, SubsetJointSvVcf keeps only that sample
    # and GATKSVPipelineSingleSample is skipped. If sample_id is absent, fall
    # back to GATKSVPipelineSingleSample. Panel / evidence inputs below are
    # unused on the prior/joint paths but remain required File inputs for the
    # single-sample path.
    String sv_joint_vcf = ""
    String sv_joint_vcf_index = ""
    # Leave unset. WDL 1.0 has no None; blank String paths coerce to this File?
    # when falling through to GATK-SV case evidence inputs.
    File? optional_file_none
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
    # When true, GatherSampleEvidence writes sample_metrics_files (manta/wham/SR/PE/counts).
    Boolean sv_run_sampleevidence_metrics = true
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

  # Blank Terra sample cells arrive as "". WDL 1.0 has no None, so use an empty
  # Array[File] when unset; Assess localizes only when this has one path.
  Array[File] phenotype_report_files = if phenotype_report != "" then [phenotype_report] else []
  # Prior File arrays are built after ProbePriorFile so a placeholder gs:// path
  # is not coerced to File (no size() 404). SV: prior → joint → GATK-SV.
  Array[File] sv_joint_vcf_files = if sv_joint_vcf != "" then [sv_joint_vcf] else []
  Array[File] sv_joint_vcf_index_files = if sv_joint_vcf_index != "" then [sv_joint_vcf_index] else []
  File? case_manta_vcf = if sv_case_manta_vcf != "" then sv_case_manta_vcf else optional_file_none
  File? case_wham_vcf = if sv_case_wham_vcf != "" then sv_case_wham_vcf else optional_file_none
  File? case_scramble_vcf = if sv_case_scramble_vcf != "" then sv_case_scramble_vcf else optional_file_none
  File? case_counts_file = if sv_case_counts_file != "" then sv_case_counts_file else optional_file_none
  File? case_pesr_disc = if sv_case_pesr_disc != "" then sv_case_pesr_disc else optional_file_none
  File? case_pesr_split = if sv_case_pesr_split != "" then sv_case_pesr_split else optional_file_none
  File? case_pesr_sd = if sv_case_pesr_sd != "" then sv_case_pesr_sd else optional_file_none

  call Tasks.FetchRepoAssets {
    input:
      git_commit = assets_git_commit,
      python_docker = python_docker
  }
  call Tasks.ProbePriorFile as ProbeIntromePrior {
    input:
      path = introme_prior,
      kind = "introme",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeVepPrior {
    input:
      path = vep_prior_vcf,
      kind = "vcf",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeVepPriorIndex {
    input:
      path = vep_prior_vcf_index,
      kind = "any",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeSmallPrior {
    input:
      path = small_vcf_prior,
      kind = "vcf",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeSmallPriorIndex {
    input:
      path = small_vcf_prior_index,
      kind = "any",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeSvPrior {
    input:
      path = sv_prior_vcf,
      kind = "vcf",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeSvPriorIndex {
    input:
      path = sv_prior_vcf_index,
      kind = "any",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeEhPrior {
    input:
      path = expansionhunter_vcf_prior,
      kind = "vcf",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeGauchianPrior {
    input:
      path = gauchian_json_prior,
      kind = "json",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeMtdnaPrior {
    input:
      path = mtdna_vcf_prior,
      kind = "vcf",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbePanelPrior {
    input:
      path = panel_csv_prior,
      kind = "csv",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeAviPrior {
    input:
      path = avi_tsv_prior,
      kind = "tsv",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  call Tasks.ProbePriorFile as ProbeGba1Prior {
    input:
      path = gba1_pd_browser_prior,
      kind = "tsv",
      probe_prior_file_py = FetchRepoAssets.probe_prior_file_py,
      python_docker = python_docker,
      preemptible = preemptible
  }
  Boolean use_prior_introme = ProbeIntromePrior.usable
  Boolean use_prior_vep = ProbeVepPrior.usable
  Boolean use_prior_small = ProbeSmallPrior.usable
  Boolean use_prior_sv = ProbeSvPrior.usable
  Boolean use_prior_eh = ProbeEhPrior.usable
  Boolean use_prior_gauchian = ProbeGauchianPrior.usable
  Boolean use_prior_mtdna = ProbeMtdnaPrior.usable
  Boolean use_prior_panel = ProbePanelPrior.usable
  Boolean use_prior_avi = ProbeAviPrior.usable
  Boolean use_prior_gba1 = ProbeGba1Prior.usable
  Array[File] introme_prior_files = if use_prior_introme then [introme_prior] else []
  Array[File] vep_prior_vcf_files = if use_prior_vep then [vep_prior_vcf] else []
  Array[File] vep_prior_vcf_index_files = if ProbeVepPriorIndex.usable then [vep_prior_vcf_index] else []
  Array[File] small_vcf_prior_files = if use_prior_small then [small_vcf_prior] else []
  Array[File] small_vcf_prior_index_files = if ProbeSmallPriorIndex.usable then [small_vcf_prior_index] else []
  Array[File] sv_prior_vcf_files = if use_prior_sv then [sv_prior_vcf] else []
  Array[File] sv_prior_vcf_index_files = if ProbeSvPriorIndex.usable then [sv_prior_vcf_index] else []
  Array[File] expansionhunter_vcf_prior_files = if use_prior_eh then [expansionhunter_vcf_prior] else []
  Array[File] gauchian_json_prior_files = if use_prior_gauchian then [gauchian_json_prior] else []
  Array[File] mtdna_vcf_prior_files = if use_prior_mtdna then [mtdna_vcf_prior] else []
  Array[File] panel_csv_prior_files = if use_prior_panel then [panel_csv_prior] else []
  Array[File] avi_tsv_prior_files = if use_prior_avi then [avi_tsv_prior] else []
  Array[File] gba1_pd_browser_prior_files = if use_prior_gba1 then [gba1_pd_browser_prior] else []
  File sample_panel_csv = select_first(flatten([panel_csv_prior_files, [FetchRepoAssets.panel_csv]]))
  File sample_gba1_pd_browser = select_first(flatten([gba1_pd_browser_prior_files, [FetchRepoAssets.gba1_pd_browser]]))

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
  call Tasks.RescueSites {
    input:
      vcf = GenotypeHardFilter.vcf,
      vcf_index = GenotypeHardFilter.vcf_index,
      rescue_nonpass_sites_py = FetchRepoAssets.rescue_nonpass_sites_py,
      clinvar_plp_vcf = clinvar_plp_vcf,
      clinvar_plp_vcf_index = clinvar_plp_vcf_index,
      gba1_pd_browser = sample_gba1_pd_browser,
      sample_id = sample_id,
      python_docker = python_docker,
      preemptible = preemptible
  }
  if (!use_prior_small) {
    call Tasks.NormPass {
      input:
        vcf = RescueSites.rescued_vcf,
        vcf_index = RescueSites.rescued_vcf_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        sample_id = sample_id,
        preemptible = preemptible
    }
  }
  if (use_prior_small) {
    call Tasks.PreferOrBuildTabix as PriorSmallIndex {
      input:
        vcf = small_vcf_prior_files[0],
        index_files = small_vcf_prior_index_files,
        name = sample_id + ".small.vcf.gz.tbi",
        preemptible = preemptible
    }
  }
  File sample_small_vcf = select_first(flatten([small_vcf_prior_files, select_all([NormPass.normalized_vcf])]))
  File sample_small_vcf_index = select_first(flatten([select_all([PriorSmallIndex.out]), select_all([NormPass.normalized_vcf_index])]))
  if (!use_prior_avi) {
    call Tasks.AnnotateAvi {
      input:
        vcf = sample_small_vcf,
        annotate_avi_py = FetchRepoAssets.annotate_avi_py,
        avi_table = avi_table,
        avi_table_index = avi_table_index,
        preemptible = preemptible
    }
  }
  Array[File] avi_tsv_files = flatten([avi_tsv_prior_files, select_all([AnnotateAvi.avi_tsv])])
  File? sample_avi_tsv = if length(avi_tsv_files) > 0 then avi_tsv_files[0] else optional_file_none
  if (!use_prior_vep) {
    call Tasks.VepAnnotate {
      input:
        vcf = sample_small_vcf,
        vcf_index = sample_small_vcf_index,
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
    if (defined(gnomad_af_vcf) && defined(gnomad_af_vcf_index)) {
      call Tasks.AnnotateGnomadAf {
        input:
          vcf = VepAnnotate.vep_vcf,
          vcf_index = VepAnnotate.vep_index,
          sites_vcf = select_first([gnomad_af_vcf]),
          sites_vcf_index = select_first([gnomad_af_vcf_index]),
          annotate_gnomad_af_py = FetchRepoAssets.annotate_gnomad_af_py,
          sample_id = sample_id,
          preemptible = preemptible
      }
    }
    call Tasks.LiftGnomadInfo {
      input:
        vcf = select_first([AnnotateGnomadAf.annotated_vcf, VepAnnotate.vep_vcf]),
        vcf_index = select_first([AnnotateGnomadAf.annotated_vcf_index, VepAnnotate.vep_index]),
        lift_gnomad_from_csq_py = FetchRepoAssets.lift_gnomad_from_csq_py,
        sample_id = sample_id,
        preemptible = preemptible
    }
    if (defined(clinvar_plp_vcf) && defined(clinvar_plp_vcf_index)) {
      call Tasks.AnnotateClinvarPlp as AnnotateClinvarSmall {
        input:
          vcf = LiftGnomadInfo.lifted_vcf,
          vcf_index = LiftGnomadInfo.lifted_vcf_index,
          sites_vcf = select_first([clinvar_plp_vcf]),
          sites_vcf_index = select_first([clinvar_plp_vcf_index]),
          annotate_clinvar_plp_py = FetchRepoAssets.annotate_clinvar_plp_py,
          sample_id = sample_id,
          python_docker = python_docker,
          preemptible = preemptible
      }
    }
  }
  if (use_prior_vep) {
    call Tasks.PreferOrBuildTabix as PriorVepIndex {
      input:
        vcf = vep_prior_vcf_files[0],
        index_files = vep_prior_vcf_index_files,
        name = sample_id + ".vep.vcf.gz.tbi",
        preemptible = preemptible
    }
  }
  File sample_vep_vcf = select_first(flatten([vep_prior_vcf_files, select_all([AnnotateClinvarSmall.annotated_vcf]), select_all([LiftGnomadInfo.lifted_vcf])]))
  File sample_vep_vcf_index = select_first(flatten([select_all([PriorVepIndex.out]), select_all([AnnotateClinvarSmall.annotated_vcf_index]), select_all([LiftGnomadInfo.lifted_vcf_index])]))
  if (!use_prior_eh) {
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
  }
  File sample_expansionhunter_vcf = select_first(flatten([expansionhunter_vcf_prior_files, select_all([ExpansionHunter.vcf])]))
  if (!use_prior_gauchian) {
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
  }
  File sample_gauchian_json = select_first(flatten([gauchian_json_prior_files, select_all([Gauchian.json])]))
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
  if (!use_prior_mtdna) {
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
    call Tasks.AnnotateMtdna {
      input:
        vcf = Mito.out_vcf,
        annotate_mtdna_py = FetchRepoAssets.annotate_mtdna_py,
        mitimpact_csv = mtdna_mitimpact,
        sample_id = sample_id,
        python_docker = python_docker,
        preemptible = preemptible
    }
  }
  File sample_mtdna_vcf = select_first(flatten([mtdna_vcf_prior_files, select_all([AnnotateMtdna.mtdna_vcf])]))
  String sample_mtdna_annotation = if use_prior_mtdna then "prior" else select_first([AnnotateMtdna.annotation_status])
  # SV precedence: prior sample VCF → joint subset → GATK-SV call (with optional
  # case_* evidence). Blank Terra strings are ignored so the next path runs.
  # scatter+if keeps Array[File] outputs for flatten with the GatkSv if-branch.
  if (use_prior_sv) {
    call Tasks.PreferOrBuildTabix as PriorSvIndex {
      input:
        vcf = sv_prior_vcf_files[0],
        index_files = sv_prior_vcf_index_files,
        name = sample_id + ".sv.vcf.gz.tbi",
        preemptible = preemptible
    }
  }
  Boolean try_joint_sv = length(sv_joint_vcf_files) > 0 && !use_prior_sv
  scatter (unused in if try_joint_sv then [1] else []) {
    call Tasks.SubsetJointSvVcf {
      input:
        joint_vcf = sv_joint_vcf_files[0],
        joint_vcf_index_files = sv_joint_vcf_index_files,
        sample_id = sample_id,
        preemptible = preemptible
    }
  }
  Boolean use_joint_sv = !use_prior_sv && select_first(flatten([SubsetJointSvVcf.sample_found, [false]]))
  # if (not scatter) so evidence intermediates are singular File? for Terra columns.
  if (!use_prior_sv && !use_joint_sv) {
    call GatkSv.GATKSVPipelineSingleSample as GatkSvSingle {
      input:
        sample_id = sample_id,
        batch = "pd-sample",
        bam_or_cram_file = ReleaseBam.released,
        bam_or_cram_index = ReleaseBai.released,
        reference_version = "38",
        use_melt = sv_use_melt,
        use_scramble = sv_use_scramble,
        case_manta_vcf = case_manta_vcf,
        case_wham_vcf = case_wham_vcf,
        case_scramble_vcf = case_scramble_vcf,
        case_counts_file = case_counts_file,
        case_pe_file = case_pesr_disc,
        case_sr_file = case_pesr_split,
        case_sd_file = case_pesr_sd,
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
        run_sampleevidence_metrics = sv_run_sampleevidence_metrics,
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
    # Ensure caller VCF indexes exist for Terra write-back even when case_* skipped GSE.
    call Tasks.PreferOrBuildTabix as EvidenceMantaIndex {
      input:
        vcf = GatkSvSingle.manta_vcf,
        index_files = select_all([GatkSvSingle.manta_index]),
        name = sample_id + ".manta.vcf.gz.tbi",
        preemptible = preemptible
    }
    call Tasks.PreferOrBuildTabix as EvidenceWhamIndex {
      input:
        vcf = GatkSvSingle.wham_vcf,
        index_files = select_all([GatkSvSingle.wham_index]),
        name = sample_id + ".wham.vcf.gz.tbi",
        preemptible = preemptible
    }
    call Tasks.PreferOrBuildTabix as EvidenceScrambleIndex {
      input:
        vcf = GatkSvSingle.scramble_vcf,
        index_files = select_all([GatkSvSingle.scramble_index]),
        name = sample_id + ".scramble.vcf.gz.tbi",
        preemptible = preemptible
    }
  }
  # Prefer prior, then GatkSv (sample absent from joint → fallback), then joint.
  # select_all turns optional if-branch File? into a 0-or-1 Array for flatten.
  Array[File] prior_sv_vcfs = sv_prior_vcf_files
  File raw_sv_vcf = select_first(flatten([prior_sv_vcfs, select_all([GatkSvSingle.final_vcf]), SubsetJointSvVcf.vcf]))
  File raw_sv_vcf_index = select_first(flatten([select_all([PriorSvIndex.out]), select_all([GatkSvSingle.final_vcf_idx]), SubsetJointSvVcf.vcf_index]))
  if (defined(clinvar_plp_vcf) && defined(clinvar_plp_vcf_index)) {
    call Tasks.AnnotateClinvarPlp as AnnotateClinvarSv {
      input:
        vcf = raw_sv_vcf,
        vcf_index = raw_sv_vcf_index,
        sites_vcf = select_first([clinvar_plp_vcf]),
        sites_vcf_index = select_first([clinvar_plp_vcf_index]),
        annotate_clinvar_plp_py = FetchRepoAssets.annotate_clinvar_plp_py,
        sample_id = sample_id + ".sv",
        python_docker = python_docker,
        preemptible = preemptible
    }
  }
  File sample_sv_vcf = select_first([AnnotateClinvarSv.annotated_vcf, raw_sv_vcf])
  File sample_sv_vcf_index = select_first([AnnotateClinvarSv.annotated_vcf_index, raw_sv_vcf_index])
  # Stage as sample_id.sv.vcf.gz (GATK-SV finals are opaque; prior/joint already match).
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
  if (!use_prior_introme) {
  call IntromeWf.Introme as Introme {
    input:
      vcf = sample_vep_vcf,
      vcf_index = sample_vep_vcf_index,
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
  }
  File sample_introme = select_first(flatten([introme_prior_files, select_all([Introme.introme_tsv])]))
  call Tasks.ReleaseAfterQc as NameIntrome {
    input:
      payload = sample_introme,
      gate = QcGate.gate,
      name = sample_id + ".introme.predictions.tsv"
  }
  call Tasks.Assess {
    input:
      sample_id = sample_id,
      sex = QcGate.sex,
      vep_vcf = sample_vep_vcf,
      small_vcf = sample_small_vcf,
      sv_vcf = NameSvVcf.released,
      sv_status = if use_prior_sv then "prior" else if use_joint_sv then "joint" else "called",
      expansionhunter_vcf = sample_expansionhunter_vcf,
      gauchian_json = sample_gauchian_json,
      mtdna_vcf = sample_mtdna_vcf,
      mtdna_annotation = sample_mtdna_annotation,
      introme_tsv = NameIntrome.released,
      introme_status = if use_prior_introme then "prior" else "called",
      avi_tsv = sample_avi_tsv,
      avi_requested = defined(avi_table) || use_prior_avi,
      qc_json = QcGate.qc_json,
      config_yaml = FetchRepoAssets.config_yaml,
      columns_yaml = FetchRepoAssets.columns_yaml,
      phenotypes_yaml = FetchRepoAssets.phenotypes_yaml,
      panel_csv = sample_panel_csv,
      gene_bed = FetchRepoAssets.gene_bed,
      exon_bed = FetchRepoAssets.exon_bed,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gtf = introme_gtf,
      gba1_pd_browser = sample_gba1_pd_browser,
      prs_weights = prs_weights,
      phenotype_report_files = phenotype_report_files,
      pdvar_docker = pdvar_docker,
      pdvar_git_url = pdvar_git_url,
      pdvar_git_commit = pdvar_git_commit
  }

  output {
    Array[File] candidates = Assess.candidates
    File vep_vcf = sample_vep_vcf
    File vep_vcf_index = sample_vep_vcf_index
    File small_vcf = sample_small_vcf
    File small_vcf_index = sample_small_vcf_index
    File sv_vcf = NameSvVcf.released
    File sv_vcf_index = NameSvVcfIndex.released
    File expansionhunter_vcf = sample_expansionhunter_vcf
    File gauchian_json = sample_gauchian_json
    File mtdna_vcf = sample_mtdna_vcf
    File introme = NameIntrome.released
    File panel_csv = sample_panel_csv
    File? avi_tsv = sample_avi_tsv
    File gba1_pd_browser = sample_gba1_pd_browser
    # GatherSampleEvidence intermediates for Terra write-back / reuse as sv_case_*.
    # File? is unset when the joint-SV path skipped GatkSvSingle. Metrics stay
    # Array[File] (empty on the joint path).
    File? coverage_counts = GatkSvSingle.coverage_counts
    File? manta_vcf = GatkSvSingle.manta_vcf
    File? manta_index = EvidenceMantaIndex.out
    File? pesr_disc = GatkSvSingle.pesr_disc
    File? pesr_disc_index = GatkSvSingle.pesr_disc_index
    File? pesr_sd = GatkSvSingle.pesr_sd
    File? pesr_sd_index = GatkSvSingle.pesr_sd_index
    File? pesr_split = GatkSvSingle.pesr_split
    File? pesr_split_index = GatkSvSingle.pesr_split_index
    File? scramble_clusters = GatkSvSingle.scramble_clusters
    File? scramble_index = EvidenceScrambleIndex.out
    File? scramble_table = GatkSvSingle.scramble_table
    File? scramble_vcf = GatkSvSingle.scramble_vcf
    File? wham_index = EvidenceWhamIndex.out
    File? wham_vcf = GatkSvSingle.wham_vcf
    Array[File] sample_metrics_files = select_first([GatkSvSingle.sample_metrics_files, []])
  }
}
