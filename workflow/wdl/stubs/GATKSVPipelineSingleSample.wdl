version 1.0

struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
}

workflow GATKSVPipelineSingleSample {
  meta {
    description: "Offline fallback only. The committed PdSample.wdl imports a vendored gatk-sv v1.1.1 GATKSVPipelineSingleSample with absolute nested imports. validate_wdl.sh --offline substitutes this file. The default path validates the real vendor import."
  }
  input {
    String sample_id
    String batch
    File bam_or_cram_file
    File bam_or_cram_index
    String? reference_version
    Boolean use_melt = false
    Boolean use_scramble = true
    Boolean use_manta = true
    Boolean use_wham = true
    File? case_manta_vcf
    File? case_wham_vcf
    File? case_scramble_vcf
    File? case_counts_file
    File? case_pe_file
    File? case_sr_file
    File? case_sd_file
    File HERVK_reference
    File LINE1_reference
    Int refine_complex_n_per_split = 15000
    File? svtk_to_gatk_script
    File ref_samples_list
    File ref_ped_file
    File genome_file
    File primary_contigs_list
    File primary_contigs_fai
    File reference_fasta
    File reference_index
    File reference_dict
    File ref_panel_vcf
    File autosome_file
    File allosome_file
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker
    String linux_docker
    String cnmops_docker
    String gatk_docker
    String condense_counts_docker
    String genomes_in_the_cloud_docker
    String samtools_cloud_docker
    String cloud_sdk_docker
    String gq_recalibrator_gatk_docker
    File preprocessed_intervals
    File manta_region_bed
    File manta_region_bed_index
    Float? manta_jobs_per_cpu
    RuntimeAttr? runtime_attr_manta
    RuntimeAttr? runtime_attr_wham
    File sd_locs_vcf
    File wham_include_list_bed_file
    Boolean run_vcf_qc
    Boolean? run_sampleevidence_metrics = false
    File wgd_scoring_mask
    Int min_svsize
    File contig_ploidy_model_tar
    File gcnv_model_tars_list
    File ref_panel_bincov_matrix
    File ref_pesr_disc_files_list
    File ref_pesr_split_files_list
    File ref_pesr_sd_files_list
    Int ref_copy_number_autosomal_contigs
    Int gcnv_qs_cutoff
    File cnmops_exclude_list
    Int matrix_qc_distance
    File ref_std_manta_vcf_tar
    File? ref_std_scramble_vcf_tar
    File ref_std_wham_vcf_tar
    File ref_panel_del_bed
    File ref_panel_dup_bed
    File depth_exclude_list
    Float depth_exclude_overlap_fraction
    Float depth_interval_overlap
    File pesr_exclude_intervals
    Float pesr_interval_overlap
    Int pesr_breakend_window
    File rmsk
    File segdups
    Int genotyping_n_per_split
    Int n_RD_genotype_bins
    File cutoffs
    File genotype_pesr_pesr_sepcutoff
    File genotype_pesr_depth_sepcutoff
    File genotype_depth_pesr_sepcutoff
    File genotype_depth_depth_sepcutoff
    File SR_metrics
    File PE_metrics
    File bin_exclude
    Float clean_vcf_min_sr_background_fail_batches
    File cytobands
    File mei_bed
    Int max_shard_size_resolve
    Int clean_vcf_max_shards_per_chrom_clean_vcf_step1
    Int clean_vcf_min_records_per_shard_clean_vcf_step1
    Int clean_vcf_samples_per_clean_vcf_step2_shard
    Int clean_vcf5_records_per_shard
    Int clean_vcf1b_records_per_shard
    File protein_coding_gtf
    File noncoding_bed
    Int annotation_sv_per_shard
    File qc_definitions
    File clustering_config_part1
    File clustering_config_part2
    File stratification_config_part1
    File stratification_config_part2
    Array[File] clustering_track_bed_files
    Array[String] clustering_track_names
    File gq_recalibrator_model_file
    Array[File] genome_tracks = []
    Array[String] recalibrate_gq_args = []
    String sl_filter_args
    String manta_docker
    String scramble_docker
    String wham_docker
  }
  output {
    File final_vcf = bam_or_cram_file
    File final_vcf_idx = bam_or_cram_index
    File coverage_counts = bam_or_cram_file
    File pesr_disc = bam_or_cram_file
    File pesr_split = bam_or_cram_file
    File pesr_sd = bam_or_cram_file
    File manta_vcf = bam_or_cram_file
    File wham_vcf = bam_or_cram_file
    File scramble_vcf = bam_or_cram_file
  }
}
