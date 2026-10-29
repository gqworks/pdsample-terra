# GATK-SV v1.1.1 inputs and Terra workspace keys

Public paths are the files `gatk-sv` tag `v1.1.1` resolves for `GATKSVPipelineSingleSample` when `scripts/inputs/build_inputs.py` is run the way Broad's single-sample test template does:

```bash
python3 scripts/inputs/build_inputs.py inputs/values inputs/templates/test/GATKSVPipelineSingleSample out \
  -a '{"single_sample":"test_single_sample_NA12878","ref_panel":"ref_panel_1kg"}'
```

`inputs/values/ref_panel_1kg.json` at v1.1.1 points at the **`gatk-sv-ref-panel-1kg-v1-0-1`** panel bundle. `inputs/values/resources_hg38.json` supplies hg38 resources (clustering/stratification configs, GQ model, genome tracks). The Terra template (`inputs/templates/terra_workspaces/single_sample/`) writes panel files to `ref_panel_*` workspace attributes and hg38 resources to `reference_*` / named twins.

`workflow/terra/PdSample.example.json` stores these `gs://` paths so `womtool validate -i` can type-check them. `workflow/terra/PdSample.terra_inputs.json` is the Terra method configuration: sample fields are `this.*`, and each GATK-SV file is `workspace.<key>`. Docker inputs are omitted so the v1.1.1 pins in `PdSample.wdl` apply.

PdSample imports a **vendored** copy of `GATKSVPipelineSingleSample.wdl` (`workflow/wdl/vendor/gatk-sv/`) that rewrites nested imports to absolute v1.1.1 URLs and pass-throughs `HERVK_reference`, `LINE1_reference`, and `RefineComplexVariants.n_per_split` so Cromwell accepts the workflow as a sub-workflow.

To skip that sub-workflow for samples in a cohort callset, set workspace attributes `sv_joint_vcf` / optional `sv_joint_vcf_index` (`workspace.sv_joint_vcf`). `SubsetJointSvVcf` keeps only `sample_id` when it is a header sample column; if the sample is absent, `GATKSVPipelineSingleSample` still runs. Panel workspace attributes above are still required by the method config for the single-sample path.

When `GATKSVPipelineSingleSample` runs, PdSample also emits GatherSampleEvidence intermediates as `Array[File]` outputs (`coverage_counts`, `manta_vcf`, `manta_index`, `pesr_disc`, `pesr_disc_index`, `pesr_sd`, `pesr_sd_index`, `pesr_split`, `pesr_split_index`, `scramble_clusters`, `scramble_index`, `scramble_table`, `scramble_vcf`, `wham_index`, `wham_vcf`, `sample_metrics_files`). Write them back to the sample table and bind the matching `sv_case_*` inputs on the next run to skip regenerating those products.

## Matched panel

Use one panel bundle. Do not mix the old `gs://gatk-sv-ref-panel-1kg/...` v0.29-beta paths with the v1.1.1 `gatk-sv-ref-panel-1kg-v1-0-1` files. Swap every `ref_panel_*` attribute together when upgrading a workspace.

Defaults: **Scramble on**, **MELT off**. Pass `ref_panel_std_scramble_vcf_tar`. Do not wire MELT unless you have a licensed melt docker and panel melt calls.

## Files (panel + new v1.1.1)

| PdSample input | v1.1.1 GATK-SV input | Public path (v1.1.1 `build_inputs`) | Terra workspace key |
| --- | --- | --- | --- |
| `sv_ref_samples_list` | `ref_samples_list` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/lists/ref_panel_1kg.samples.list` | `workspace.ref_panel_samples_list` |
| `sv_ref_ped_file` | `ref_ped_file` | `gs://gcp-public-data--broad-references/hg38/v0/sv-resources/ref-panel/1KG/v1/ped/1kg_ref_panel_v1.ped` | `workspace.ref_panel_ped_file` |
| `sv_ref_panel_vcf` | `ref_panel_vcf` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/6f760c00-1c94-4b9f-9785-43533af3a3ed/CleanVcf/.../all_batches.cleaned.vcf.gz` | `workspace.ref_panel_vcf` |
| `sv_ref_panel_bincov_matrix` | `ref_panel_bincov_matrix` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/fc210106-.../all_samples.RD.txt.gz` | `workspace.ref_panel_bincov_matrix` |
| `sv_ref_pesr_disc_files_list` | `ref_pesr_disc_files_list` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/lists/ref_panel_1kg.PE_files.list` | `workspace.ref_panel_PE_files_list` |
| `sv_ref_pesr_split_files_list` | `ref_pesr_split_files_list` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/lists/ref_panel_1kg.SR_files.list` | `workspace.ref_panel_SR_files_list` |
| `sv_ref_pesr_sd_files_list` | `ref_pesr_sd_files_list` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/lists/ref_panel_1kg.SD_files.list` | `workspace.ref_panel_SD_files_list` |
| `sv_contig_ploidy_model_tar` | `contig_ploidy_model_tar` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/b4afacf6-.../all_samples-contig-ploidy-model.tar.gz` | `workspace.ref_panel_contig_ploidy_model_tar` |
| `sv_gcnv_model_tars_list` | `gcnv_model_tars_list` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/lists/ref_panel_1kg.gcnv_model_tars.list` | `workspace.ref_panel_gcnv_model_tars_list` |
| `sv_ref_std_manta_vcf_tar` | `ref_std_manta_vcf_tar` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/fc210106-.../all_samples.manta_std.tar.gz` | `workspace.ref_panel_std_manta_vcf_tar` |
| `sv_ref_std_scramble_vcf_tar` | `ref_std_scramble_vcf_tar` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/fc210106-.../all_samples.scramble_std.tar.gz` | `workspace.ref_panel_std_scramble_vcf_tar` |
| `sv_ref_std_wham_vcf_tar` | `ref_std_wham_vcf_tar` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/fc210106-.../all_samples.wham_std.tar.gz` | `workspace.ref_panel_std_wham_vcf_tar` |
| `sv_ref_panel_del_bed` | `ref_panel_del_bed` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/fc210106-.../all_samples.DEL.bed.gz` | `workspace.ref_panel_del_bed` |
| `sv_ref_panel_dup_bed` | `ref_panel_dup_bed` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/fc210106-.../all_samples.DUP.bed.gz` | `workspace.ref_panel_dup_bed` |
| `sv_cutoffs` | `cutoffs` | `gs://gatk-sv-ref-panel-1kg-v1-0-1/submissions/9f2af2ea-.../all_samples.cutoffs` | `workspace.ref_panel_cutoffs` |
| `sv_genotype_*_sepcutoff` / `sv_*_metrics` | same names | under `gatk-sv-ref-panel-1kg-v1-0-1/submissions/6b2ae6b8-...` | matching `workspace.ref_panel_*` |
| `sv_qc_definitions` | `qc_definitions` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/single_sample.qc_definitions.tsv` | `workspace.ref_panel_qc_definitions` |
| `sv_clustering_config_part1` | `clustering_config_part1` | `gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/clustering_config.part_one.tsv` | `workspace.clustering_config_part1` |
| `sv_clustering_config_part2` | `clustering_config_part2` | `.../clustering_config.part_two.tsv` | `workspace.clustering_config_part2` |
| `sv_stratification_config_part1` | `stratification_config_part1` | `.../stratify_config.part_one.tsv` | `workspace.stratification_config_part1` |
| `sv_stratification_config_part2` | `stratification_config_part2` | `.../stratify_config.part_two.tsv` | `workspace.stratification_config_part2` |
| `sv_gq_recalibrator_model_file` | `gq_recalibrator_model_file` | `.../gatk-sv-recalibrator.hgsvc_release3.v1.model` | `workspace.recalibrate_gq_model_file` |
| `sv_mei_bed` | `mei_bed` | `.../hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz` | `workspace.reference_mei_bed` |
| `sv_protein_coding_gtf` | `protein_coding_gtf` | `.../gencode.v47.basic.protein_coding.canonical.gtf` | `workspace.reference_protein_coding_gtf` |

Full exact URIs are in `PdSample.example.json`. hg38 reference twins (`sv_genome_file`, FASTA, contigs, exclude lists, etc.) keep the same Broad public paths as before; see the example JSON.

Removed vs v0.29-beta: `sv_empty_file` / `empty_file` (not an input at v1.1.1).

## Caller flags

| PdSample input | GATK-SV input | Default |
| --- | --- | --- |
| `sv_use_melt` | `use_melt` | `false` |
| `sv_use_scramble` | `use_scramble` | `true` |

## Dockers

Leave unset in Terra so `PdSample.wdl` keeps the v1.1.1 `inputs/values/dockers.json` pins.

| PdSample input | Image (v1.1.1) |
| --- | --- |
| `sv_manta_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/manta:2023-09-14-v0.28.3-beta-3f22f94d` |
| `sv_scramble_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/scramble:2024-10-25-v0.29-beta-5ea22a52` |
| `sv_wham_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/wham:2024-10-25-v0.29-beta-5ea22a52` |
| `sv_sv_base_mini_docker` / `sv_sv_base_docker` | `.../sv-base-mini` / `sv-base` `:2024-10-25-v0.29-beta-5ea22a52` |
| `sv_sv_pipeline_docker` / `sv_sv_pipeline_qc_docker` | `.../sv-pipeline:2025-10-02-v1.1-483973d6` |
| `sv_cnmops_docker` | `.../cnmops:2025-09-02-v1.0.5-f091af0b` |
| `sv_gatk_docker` | `.../gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT` |
| `sv_gq_recalibrator_gatk_docker` | `us.gcr.io/broad-dsde-methods/markw/gatk:mw-tb-form-sv-filter-training-data-eb62a3` |
| `sv_samtools_cloud_docker` | `.../samtools-cloud:2024-10-25-v0.29-beta-5ea22a52` |

Array defaults (`sv_clustering_track_*`, `sv_genome_tracks`, `sv_recalibrate_gq_args`) and `sv_sl_filter_args` live in `PdSample.wdl` with the Broad test-template values.
