# GATK-SV v0.29-beta inputs and Terra workspace keys

Public paths are the files `gatk-sv` tag `v0.29-beta` resolves for `GATKSVPipelineSingleSample` when `scripts/inputs/build_inputs.py` is run the way `scripts/inputs/build_default_inputs.sh` does for a single sample:

```bash
python3 scripts/inputs/build_inputs.py inputs/values inputs/templates/test/GATKSVPipelineSingleSample out.json \
  -a '{"single_sample":"test_single_sample_NA12878","ref_panel":"ref_panel_1kg"}'
```

`inputs/values/ref_panel_1kg.json` supplies the 1KG panel. `inputs/values/resources_hg38.json` supplies the hg38 resource files. The single-sample Terra template (`inputs/templates/terra_workspaces/single_sample/`) writes panel files to `ref_panel_*` workspace attributes and hg38 resources to the `reference_*` twins. `qc_definitions` is the exception: the template sets `ref_panel_qc_definitions` from `resources_hg38.single_sample_qc_definitions`, and that is the file the resolved workflow input uses.

`workflow/terra/PdSample.example.json` stores these `gs://` paths so `womtool validate -i` can type-check them. `workflow/terra/PdSample.terra_inputs.json` is the Terra method configuration: sample fields are `this.sample_id`, `this.cram`, `this.crai`, and `this.sex`, and each GATK-SV file is `workspace.<key>`. Scalar tuning values that are not workspace data are the literals from the v0.29-beta `GATKSVPipelineSingleSample` test template. Docker inputs are omitted from both JSON files so the v0.29-beta pins in `PdSample.wdl` apply. Do not bind the workspace `*_docker` attributes: a featured GATK-SV workspace can mix image tags (some v0.29-beta, with `sv_pipeline_docker` / `sv_pipeline_qc_docker` at other releases).

## Matched panel

Use one panel bundle. Do not pair a cohort `ref_panel_vcf` or ped with the 1KG gCNV models or bincov matrix. The public default is the whole `ref_panel_1kg` bundle from tag v0.29-beta, including these path families that the tag already combines:

- `contig_ploidy_model_tar` is the public 1KG **v2** gCNV tar (`.../ref-panel/1KG/v2/gcnv/ref_panel_1kg_v2-contig-ploidy-model.tar.gz`). The `gcnv_model_tars` array in the same JSON (285 shards) is those public v2 model tars. `gcnv_model_tars_list`, `samples_list`, the bincov matrix, del/dup beds, cutoffs, genotype separation cutoffs, and PE/SR metrics live under batch `38c65ca4-2a07-4805-86b6-214696075fef`.
- std Manta and Wham tars are `GatherBatchEvidence/b705e9f6-c401-4a9a-a059-1a68d1ca5f98`.
- PE/SR/SD file lists are `tws-no-cram-conversion/GatherSampleEvidenceBatch`.
- The cleaned panel VCF is `mw-vcf-reshard/ref_panel_1kg.cleaned.vcf.gz`. The ped is the public v1 file `1kg_ref_panel_v1.ped`.
- `ref_panel_1kg.json` `qc_definitions` is the batch TSV under `38c65ca4-...`. The single-sample workflow input is the public v2 `single_sample.qc_definitions.tsv`, which is what `workspace.ref_panel_qc_definitions` is filled with.

That v2 ploidy tar next to batch `38c65ca4` is the combination v0.29-beta ships in `ref_panel_1kg.json`. It is not a later workspace edit. Swap panels by replacing the whole set of `ref_panel_*` attributes together, and keep `reference_*` on the hg38 resource bundle those models were trained with. The panel has to be built with a GATK-SV version compatible with v0.29-beta.

To point Terra at workspace data instead of the baked-in `gs://` paths, import `PdSample.terra_inputs.json` or set each input in the Terra Inputs tab to the expression in the last column. Sample columns stay on the sample table (`this.sample_id`, `this.cram`, `this.crai`, `this.sex`).

## Files

| PdSample input | v0.29-beta GATK-SV input | Public 1KG / hg38 path | Terra workspace key |
| --- | --- | --- | --- |
| `sv_ref_samples_list` | `ref_samples_list` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/samples_list.txt` | `workspace.ref_panel_samples_list` |
| `sv_ref_ped_file` | `ref_ped_file` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/ref-panel/1KG/v1/ped/1kg_ref_panel_v1.ped` | `workspace.ref_panel_ped_file` |
| `sv_ref_panel_vcf` | `ref_panel_vcf` | `gs://gatk-sv-ref-panel-1kg/outputs/mw-vcf-reshard/ref_panel_1kg.cleaned.vcf.gz` | `workspace.ref_panel_vcf` |
| `sv_ref_panel_bincov_matrix` | `ref_panel_bincov_matrix` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-EvidenceQC/EvidenceQC/f7f8f407-0a3c-4c3d-982f-7bd181ac965a/call-MakeBincovMatrix/MakeBincovMatrix/0106cefa-5ebf-4131-9073-15413bf93a83/call-ZPaste/cacheCopy/ref_panel_1kg.RD.txt.gz` | `workspace.ref_panel_bincov_matrix` |
| `sv_ref_pesr_disc_files_list` | `ref_pesr_disc_files_list` | `gs://gatk-sv-ref-panel-1kg/outputs/tws-no-cram-conversion/GatherSampleEvidenceBatch/pe_files_list.txt` | `workspace.ref_panel_PE_files_list` |
| `sv_ref_pesr_split_files_list` | `ref_pesr_split_files_list` | `gs://gatk-sv-ref-panel-1kg/outputs/tws-no-cram-conversion/GatherSampleEvidenceBatch/sr_files_list.txt` | `workspace.ref_panel_SR_files_list` |
| `sv_ref_pesr_sd_files_list` | `ref_pesr_sd_files_list` | `gs://gatk-sv-ref-panel-1kg/outputs/tws-no-cram-conversion/GatherSampleEvidenceBatch/sd_files_list.txt` | `workspace.ref_panel_SD_files_list` |
| `sv_contig_ploidy_model_tar` | `contig_ploidy_model_tar` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/gcnv/ref_panel_1kg_v2-contig-ploidy-model.tar.gz` | `workspace.ref_panel_contig_ploidy_model_tar` |
| `sv_gcnv_model_tars_list` | `gcnv_model_tars_list` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/gcnv_model_tars_list.txt` | `workspace.ref_panel_gcnv_model_tars_list` |
| `sv_ref_std_manta_vcf_tar` | `ref_std_manta_vcf_tar` | `gs://gatk-sv-ref-panel-1kg/outputs/GatherBatchEvidence/b705e9f6-c401-4a9a-a059-1a68d1ca5f98/ref_panel_1kg.manta_std.tar.gz` | `workspace.ref_panel_std_manta_vcf_tar` |
| `sv_ref_std_wham_vcf_tar` | `ref_std_wham_vcf_tar` | `gs://gatk-sv-ref-panel-1kg/outputs/GatherBatchEvidence/b705e9f6-c401-4a9a-a059-1a68d1ca5f98/ref_panel_1kg.wham_std.tar.gz` | `workspace.ref_panel_std_wham_vcf_tar` |
| `sv_ref_panel_del_bed` | `ref_panel_del_bed` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GATKSVPipelinePhase1/GATKSVPipelinePhase1/acce2c71-7458-4205-ae13-624f6efc9956/call-GatherBatchEvidence/GatherBatchEvidence/366e817a-feb2-4d3c-915c-a8bd25529b81/call-MergeDepth/MergeDepth/b58e04ee-846d-4be6-bba3-6dda2f2c995e/call-MergeSet_del/cacheCopy/ref_panel_1kg.DEL.bed.gz` | `workspace.ref_panel_del_bed` |
| `sv_ref_panel_dup_bed` | `ref_panel_dup_bed` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GATKSVPipelinePhase1/GATKSVPipelinePhase1/acce2c71-7458-4205-ae13-624f6efc9956/call-GatherBatchEvidence/GatherBatchEvidence/366e817a-feb2-4d3c-915c-a8bd25529b81/call-MergeDepth/MergeDepth/b58e04ee-846d-4be6-bba3-6dda2f2c995e/call-MergeSet_dup/cacheCopy/ref_panel_1kg.DUP.bed.gz` | `workspace.ref_panel_dup_bed` |
| `sv_cutoffs` | `cutoffs` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GATKSVPipelinePhase1/GATKSVPipelinePhase1/acce2c71-7458-4205-ae13-624f6efc9956/call-FilterBatch/FilterBatch/184defa3-e61c-4757-9962-f685f6d0d204/call-FilterBatchSites/FilterBatchSites/13801c7d-9478-40dc-9a4a-2e80cbab1136/call-AdjudicateSV/cacheCopy/ref_panel_1kg.cutoffs` | `workspace.ref_panel_cutoffs` |
| `sv_genotype_pesr_pesr_sepcutoff` | `genotype_pesr_pesr_sepcutoff` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GenotypeBatch/GenotypeBatch/ad17f522-0950-4f0a-9148-a13f689082ed/call-GenotypePESRPart1/GenotypePESRPart1/40ec6d76-dd1c-432d-bfab-bc4426d0b1ec/call-TrainRDGenotyping/TrainRDGenotyping/e5540a96-9072-4719-bcfb-afccdfec15c6/call-UpdateCutoff/cacheCopy/ref_panel_1kg.pesr.pesr_sepcutoff.txt` | `workspace.ref_panel_genotype_pesr_pesr_sepcutoff` |
| `sv_genotype_pesr_depth_sepcutoff` | `genotype_pesr_depth_sepcutoff` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GenotypeBatch/GenotypeBatch/ad17f522-0950-4f0a-9148-a13f689082ed/call-GenotypePESRPart1/GenotypePESRPart1/40ec6d76-dd1c-432d-bfab-bc4426d0b1ec/call-TrainRDGenotyping/TrainRDGenotyping/e5540a96-9072-4719-bcfb-afccdfec15c6/call-UpdateCutoff/cacheCopy/ref_panel_1kg.pesr.depth_sepcutoff.txt` | `workspace.ref_panel_genotype_pesr_depth_sepcutoff` |
| `sv_genotype_depth_pesr_sepcutoff` | `genotype_depth_pesr_sepcutoff` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GenotypeBatch/GenotypeBatch/ad17f522-0950-4f0a-9148-a13f689082ed/call-GenotypeDepthPart1/GenotypeDepthPart1/be787bb1-25ff-4a54-92a6-2fa6faaee8ec/call-TrainRDGenotyping/TrainRDGenotyping/aa7bc7f2-9779-4242-aa5d-8a8ea4375af9/call-UpdateCutoff/cacheCopy/ref_panel_1kg.depth.pesr_sepcutoff.txt` | `workspace.ref_panel_genotype_depth_pesr_sepcutoff` |
| `sv_genotype_depth_depth_sepcutoff` | `genotype_depth_depth_sepcutoff` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GenotypeBatch/GenotypeBatch/ad17f522-0950-4f0a-9148-a13f689082ed/call-GenotypeDepthPart1/GenotypeDepthPart1/be787bb1-25ff-4a54-92a6-2fa6faaee8ec/call-TrainRDGenotyping/TrainRDGenotyping/aa7bc7f2-9779-4242-aa5d-8a8ea4375af9/call-UpdateCutoff/cacheCopy/ref_panel_1kg.depth.depth_sepcutoff.txt` | `workspace.ref_panel_genotype_depth_depth_sepcutoff` |
| `sv_SR_metrics` | `SR_metrics` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GenotypeBatch/GenotypeBatch/ad17f522-0950-4f0a-9148-a13f689082ed/call-GenotypePESRPart1/GenotypePESRPart1/40ec6d76-dd1c-432d-bfab-bc4426d0b1ec/call-TrainSRGenotyping/TrainSRGenotyping/e81a3092-0c6e-4768-afed-27d7018175dc/call-GenotypeSRPart1/cacheCopy/ref_panel_1kg.sr_metric_file.txt` | `workspace.ref_panel_SR_metrics` |
| `sv_PE_metrics` | `PE_metrics` | `gs://gatk-sv-ref-panel-1kg/outputs/GATKSVPipelineBatch/38c65ca4-2a07-4805-86b6-214696075fef/call-GenotypeBatch/GenotypeBatch/ad17f522-0950-4f0a-9148-a13f689082ed/call-GenotypePESRPart1/GenotypePESRPart1/40ec6d76-dd1c-432d-bfab-bc4426d0b1ec/call-TrainPEGenotyping/TrainPEGenotyping/8c1271fd-027d-4669-b516-b41ee77d8997/call-GenotypePEPart1/cacheCopy/ref_panel_1kg.pe_metric_file.txt` | `workspace.ref_panel_PE_metrics` |
| `sv_qc_definitions` | `qc_definitions` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/single_sample.qc_definitions.tsv` | `workspace.ref_panel_qc_definitions` |
| `sv_genome_file` | `genome_file` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/hg38.genome` | `workspace.reference_genome_file` |
| `sv_primary_contigs_list` | `primary_contigs_list` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs.list` | `workspace.reference_primary_contigs_list` |
| `sv_primary_contigs_fai` | `primary_contigs_fai` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/contig.fai` | `workspace.reference_primary_contigs_fai` |
| `sv_reference_fasta` | `reference_fasta` | `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta` | `workspace.reference_fasta` |
| `sv_reference_index` | `reference_index` | `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai` | `workspace.reference_index` |
| `sv_reference_dict` | `reference_dict` | `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict` | `workspace.reference_dict` |
| `sv_autosome_file` | `autosome_file` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/autosome.fai` | `workspace.reference_autosome_file` |
| `sv_allosome_file` | `allosome_file` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/allosome.fai` | `workspace.reference_allosome_file` |
| `sv_preprocessed_intervals` | `preprocessed_intervals` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/preprocessed_intervals.interval_list` | `workspace.reference_preprocessed_intervals` |
| `sv_manta_region_bed` | `manta_region_bed` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz` | `workspace.reference_manta_region_bed` |
| `sv_manta_region_bed_index` | `manta_region_bed_index` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz.tbi` | `workspace.reference_manta_region_bed_index` |
| `sv_sd_locs_vcf` | `sd_locs_vcf` | `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf` | `workspace.reference_sd_locs_vcf` |
| `sv_wham_include_list_bed_file` | `wham_include_list_bed_file` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/wham_whitelist.bed` | `workspace.reference_wham_include_list_bed_file` |
| `sv_wgd_scoring_mask` | `wgd_scoring_mask` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/wgd_scoring_mask.hg38.gnomad_v3.bed` | `workspace.reference_wgd_scoring_mask` |
| `sv_cnmops_exclude_list` | `cnmops_exclude_list` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/GRCh38_Nmask.bed` | `workspace.reference_cnmops_exclude_list` |
| `sv_depth_exclude_list` | `depth_exclude_list` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/depth_blacklist.sorted.bed.gz` | `workspace.reference_depth_exclude_list` |
| `sv_pesr_exclude_intervals` | `pesr_exclude_intervals` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed.gz` | `workspace.reference_pesr_exclude_list` |
| `sv_rmsk` | `rmsk` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/hg38.randomForest_blacklist.withRepMask.bed.gz` | `workspace.reference_rmsk` |
| `sv_segdups` | `segdups` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/hg38.SD_gaps_Cen_Tel_Heter_Satellite_lumpy.blacklist.sorted.merged.bed.gz` | `workspace.reference_segdups` |
| `sv_bin_exclude` | `bin_exclude` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/bin_exclude.hg38.gatkcov.bed.gz` | `workspace.reference_bin_exclude` |
| `sv_cytobands` | `cytobands` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/cytobands_hg38.bed.gz` | `workspace.reference_cytobands` |
| `sv_mei_bed` | `mei_bed` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/mei_hg38.bed.gz` | `workspace.reference_mei_bed` |
| `sv_protein_coding_gtf` | `protein_coding_gtf` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/MANE.GRCh38.v1.2.ensembl_genomic.gtf` | `workspace.reference_protein_coding_gtf` |
| `sv_noncoding_bed` | `noncoding_bed` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed` | `workspace.reference_noncoding_bed` |
| `sv_empty_file` | `empty_file` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/empty.file` | `workspace.reference_empty_file` |

## Scalars

Literals are the v0.29-beta test template values (`inputs/templates/test/GATKSVPipelineSingleSample`, resolved by `build_inputs.py`). `ref_copy_number_autosomal_contigs` is also workspace data (`reference_copy_number_autosomal_contigs`, value 2 in `resources_hg38.json`). The example JSON keeps the number `2` so womtool sees an Int. The Terra inputs file uses the workspace expression.

| PdSample input | v0.29-beta GATK-SV input | Template value | Terra expression |
| --- | --- | --- | --- |
| `sv_run_vcf_qc` | `run_vcf_qc` | `false` | `false` |
| `sv_min_svsize` | `min_svsize` | `50` | `50` |
| `sv_ref_copy_number_autosomal_contigs` | `ref_copy_number_autosomal_contigs` | `2` | `workspace.reference_copy_number_autosomal_contigs` |
| `sv_gcnv_qs_cutoff` | `gcnv_qs_cutoff` | `30` | `30` |
| `sv_matrix_qc_distance` | `matrix_qc_distance` | `1000000` | `1000000` |
| `sv_depth_exclude_overlap_fraction` | `depth_exclude_overlap_fraction` | `0.5` | `0.5` |
| `sv_depth_interval_overlap` | `depth_interval_overlap` | `0.8` | `0.8` |
| `sv_pesr_interval_overlap` | `pesr_interval_overlap` | `0.1` | `0.1` |
| `sv_pesr_breakend_window` | `pesr_breakend_window` | `300` | `300` |
| `sv_genotyping_n_per_split` | `genotyping_n_per_split` | `1000` | `1000` |
| `sv_n_RD_genotype_bins` | `n_RD_genotype_bins` | `100000` | `100000` |
| `sv_clean_vcf_min_sr_background_fail_batches` | `clean_vcf_min_sr_background_fail_batches` | `0.5` | `0.5` |
| `sv_max_shard_size_resolve` | `max_shard_size_resolve` | `500` | `500` |
| `sv_clean_vcf_max_shards_per_chrom_clean_vcf_step1` | `clean_vcf_max_shards_per_chrom_clean_vcf_step1` | `200` | `200` |
| `sv_clean_vcf_min_records_per_shard_clean_vcf_step1` | `clean_vcf_min_records_per_shard_clean_vcf_step1` | `5000` | `5000` |
| `sv_clean_vcf_samples_per_clean_vcf_step2_shard` | `clean_vcf_samples_per_clean_vcf_step2_shard` | `100` | `100` |
| `sv_clean_vcf5_records_per_shard` | `clean_vcf5_records_per_shard` | `5000` | `5000` |
| `sv_clean_vcf1b_records_per_shard` | `clean_vcf1b_records_per_shard` | `10000` | `10000` |
| `sv_annotation_sv_per_shard` | `annotation_sv_per_shard` | `5000` | `5000` |

## Dockers

Leave these unset in Terra so `PdSample.wdl` keeps the v0.29-beta `inputs/values/dockers.json` pins. Workspace keys of the same name exist and are the wrong place to bind them.

| PdSample input | v0.29-beta GATK-SV input | Pinned image |
| --- | --- | --- |
| `sv_manta_docker` | `manta_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/manta:2023-09-14-v0.28.3-beta-3f22f94d` |
| `sv_scramble_docker` | `scramble_docker` | `us.gcr.io/broad-dsde-methods/markw/scramble:mw-scramble-99af4c50` |
| `sv_wham_docker` | `wham_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/wham:2024-01-24-v0.28.4-beta-9debd6d7` |
| `sv_sv_base_mini_docker` | `sv_base_mini_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-01-24-v0.28.4-beta-9debd6d7` |
| `sv_sv_base_docker` | `sv_base_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/sv-base:2024-01-24-v0.28.4-beta-9debd6d7` |
| `sv_sv_pipeline_docker` | `sv_pipeline_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-08-19-v0.28.5-beta-84a0627d` |
| `sv_sv_pipeline_qc_docker` | `sv_pipeline_qc_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-08-19-v0.28.5-beta-84a0627d` |
| `sv_linux_docker` | `linux_docker` | `marketplace.gcr.io/google/ubuntu1804` |
| `sv_cnmops_docker` | `cnmops_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/cnmops:2024-06-04-v0.28.5-beta-a8dfecba` |
| `sv_gatk_docker` | `gatk_docker` | `us.gcr.io/broad-dsde-methods/eph/gatk:2024-07-02-4.6.0.0-1-g4af2b49e9-NIGHTLY-SNAPSHOT` |
| `sv_condense_counts_docker` | `condense_counts_docker` | `us.gcr.io/broad-dsde-methods/tsharpe/gatk:4.2.6.1-57-g9e03432` |
| `sv_genomes_in_the_cloud_docker` | `genomes_in_the_cloud_docker` | `us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135` |
| `sv_samtools_cloud_docker` | `samtools_cloud_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/samtools-cloud:2024-01-24-v0.28.4-beta-9debd6d7` |
| `sv_cloud_sdk_docker` | `cloud_sdk_docker` | `google/cloud-sdk` |

Optional v0.29-beta inputs that `PdSample` does not pass (`chr_x`, `chr_y`, `external_af_ref_bed`, `external_af_ref_bed_prefix`, `AnnotateVcf.par_bed`, `contig_ploidy_priors`, `seed_cutoffs`, `ref_std_melt_vcf_tar`, `max_ref_panel_carrier_freq`, `clean_vcf_random_seed`) keep the upstream defaults. `HERVK_reference`, `LINE1_reference`, clustering and stratification configs, and the GQ recalibrator model are not inputs of `GATKSVPipelineSingleSample` at this tag.
