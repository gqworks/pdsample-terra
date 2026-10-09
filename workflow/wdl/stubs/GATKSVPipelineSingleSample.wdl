version 1.0

workflow GATKSVPipelineSingleSample {
  meta {
    description: "Offline fallback only. The committed PdSample.wdl imports gatk-sv v0.29-beta. validate_wdl.sh --offline substitutes this file. The default path validates the real import."
  }
  input {
    String sample_id
    String batch
    File bam_or_cram_file
    File bam_or_cram_index
    String reference_version
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
    File preprocessed_intervals
    File manta_region_bed
    File manta_region_bed_index
    File sd_locs_vcf
    File wham_include_list_bed_file
    Boolean run_vcf_qc
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
    File ref_std_wham_vcf_tar
    File ref_panel_del_bed
    File ref_panel_dup_bed
    File depth_exclude_list
    File empty_file
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
    String manta_docker
    String scramble_docker
    String wham_docker
  }
  output {
    File final_vcf = bam_or_cram_file
    File final_vcf_idx = bam_or_cram_index
  }
}
