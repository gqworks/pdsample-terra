version 1.0

workflow MitochondriaPipeline {
  meta {
    description: "Validation stub for miniwdl. The committed PdSample.wdl imports GATK 4.7.0.0. womtool validates that real import. miniwdl rejects optional arithmetic in the upstream AlignAndCall.wdl."
  }
  input {
    File wgs_aligned_input_bam_or_cram
    File wgs_aligned_input_bam_or_cram_index
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict
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
    String? requester_pays_project
    Boolean compress_output_vcf = false
    Int? preemptible_tries
  }
  output {
    File out_vcf = mt_fasta
    File out_vcf_index = mt_fasta_index
  }
}
