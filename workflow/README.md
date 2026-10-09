# PdSample workflow

WDL 1.0 for Terra / Cromwell. Entry workflow: `workflow/wdl/PdSample.wdl`.

Validate every WDL with:

```bash
bash workflow/scripts/validate_wdl.sh
```

The default path runs `womtool validate` on the committed files, so both remote imports are resolved for real:

`https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.29-beta/wdl/GATKSVPipelineSingleSample.wdl`

and

`https://raw.githubusercontent.com/broadinstitute/gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl`.

`miniwdl check` uses that same GATK-SV import. It substitutes `stubs/MitochondriaPipeline.wdl` because upstream `AlignAndCall.wdl` multiplies an optional `Int?` and miniwdl rejects that. `womtool validate` of the real mitochondria import succeeds, including as a sub-workflow. `nvidiaDriverVersion` is a Cromwell runtime attribute; miniwdl warns that it does not know the key. That warning is expected on the Introme GPU path and does not fail the check.

`bash workflow/scripts/validate_wdl.sh --offline` is the only path that replaces the GATK-SV import with `stubs/GATKSVPipelineSingleSample.wdl`. Use it when GitHub is unreachable. It is not the Dockstore or Terra check.

## GATK-SV tag

`GATKSVPipelineSingleSample` is pinned to **v0.29-beta**, not v1.1.1. v1.1.1 (and every release from v1.0 through v1.1.1, plus `main` as of this pin) calls `MakeCohortVcf` and `RefineComplexVariants` without declaring or passing `HERVK_reference`, `LINE1_reference`, and `n_per_split`. Those nested inputs are required and have no defaults, so Cromwell rejects the single-sample workflow as a sub-workflow:

`To be called as a sub-workflow it must declare and pass-through the following values via workflow inputs: MakeCohortVcf.HERVK_reference, MakeCohortVcf.LINE1_reference, RefineComplexVariants.n_per_split`

v0.29-beta is the newest release tag where those inputs are not required inside the nested workflows, so `womtool validate` accepts `PdSample.wdl` with the real import. Docker defaults are the images in that tag's `inputs/values/dockers.json`. v0.29-beta also requires `empty_file` (passed from `sv_empty_file`) and does not take the v1.1.1 clustering-config, GQ-recalibrator, or `sl_filter_args` inputs.

GATK-SV is a required call. Every file input of `GATKSVPipelineSingleSample` at this tag is a required `PdSample` input. `workflow/terra/PdSample.example.json` fills them from the public 1KG panel and hg38 resources that tag resolves. `workflow/terra/PdSample.terra_inputs.json` is the same method config with sample columns as `this.*` and panel/resource files as `workspace.<key>`. The key table is `workflow/terra/gatksv_panel_keys.md`.

Leave the `sv_*_docker` inputs unset. A featured GATK-SV Terra workspace stores docker attributes at mixed versions (some v0.29-beta; `sv_pipeline_docker` and `sv_pipeline_qc_docker` can be an older or newer tag). `PdSample` keeps the v0.29-beta pins instead of those attributes.

The public `ref_panel_1kg` bundle already mixes a 1KG v2 contig-ploidy model tar with batch `38c65ca4-2a07-4805-86b6-214696075fef` for the sample list, bincov matrix, beds, and cutoffs. That is what `inputs/values/ref_panel_1kg.json` ships. Use the whole bundle together. Details are in the key table.

## Task graph

| Step | Image | Version |
| --- | --- | --- |
| CRAM header | `quay.io/biocontainers/samtools:1.24--h9dcdb79_1` | samtools 1.24 |
| Reference dict, HaplotypeCaller, hard filters | `us.gcr.io/broad-gatk/gatk:4.6.1.0` | GATK 4.6.1.0 |
| Reference check, QC summary | `python:3.11.11-slim-bookworm` | stdlib scripts |
| VerifyBamID2 | `quay.io/biocontainers/verifybamid2:2.0.3--hc004090_0` | 2.0.3 |
| somalier sex | `quay.io/biocontainers/somalier:0.3.5--h5205c93_0` | 0.3.5 |
| mosdepth | `quay.io/biocontainers/mosdepth:0.3.14--h87be163_2` | 0.3.14 |
| PASS and norm | `quay.io/biocontainers/bcftools:1.24--h118bc1c_4` | bcftools 1.24 |
| VEP | `quay.io/biocontainers/ensembl-vep:116.2--pl5321h2a3209d_0` | 116.2 |
| ExpansionHunter | `quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5` | 5.0.0, 19-locus catalog |
| Gauchian | `python:3.11.11-slim-bookworm` unless `gauchian_docker` is set | pip `gauchian==1.0.2` |
| mtDNA | images inside MitochondriaPipeline | GATK 4.7.0.0 |
| GATK-SV | pins in `PdSample.wdl` (`sv_*_docker`) | gatk-sv v0.29-beta, required |
| Introme preprocessing | `gabyou/data_prepocessing:v2.0` | upstream tag, including the image-name spelling |
| VariantInfo | `gabyou/variant_info:v3.2` | v3.2 |
| SpliceAI | `headoncollusion/spliceai:v1.4-cpu` or `v1.4-gpu` | v1.4 |
| MMSplice | `headoncollusion/mmsplice:v2.2-cpu` or `v2.2-gpu` | v2.2 |
| Pangolin | `headoncollusion/pangolin:v1.3-cpu` or `v1.3-gpu` | v1.3 |
| SPiP | `headoncollusion/spip:v1.1` | v1.1 |
| Spliceogen | `headoncollusion/spliceogen:v2.1` | v2.1 |
| Introme functions / anno / infer | `headoncollusion/ag_check:v1.6` | v1.6 |
| Introme asset fetch | `python:3.11.11-slim-bookworm` | archive `d6da48bc` |
| Assess | `python:3.11.11-slim-bookworm` or your `pdvar` image | git commit or the Dockerfile |

GATK-SV image pins (overridable), from v0.29-beta `inputs/values/dockers.json`:

- manta `us.gcr.io/broad-dsde-methods/gatk-sv/manta:2023-09-14-v0.28.3-beta-3f22f94d`
- scramble `us.gcr.io/broad-dsde-methods/markw/scramble:mw-scramble-99af4c50`
- wham / sv-base / sv-base-mini / samtools-cloud `2024-01-24-v0.28.4-beta-9debd6d7`
- sv-pipeline `2024-08-19-v0.28.5-beta-84a0627d`
- cnmops `2024-06-04-v0.28.5-beta-a8dfecba`
- gatk `us.gcr.io/broad-dsde-methods/eph/gatk:2024-07-02-4.6.0.0-1-g4af2b49e9-NIGHTLY-SNAPSHOT`
- cloud SDK `google/cloud-sdk`

## Assess image

```bash
docker build -f workflow/containers/pdvar/Dockerfile -t REGISTRY/pdvar:COMMIT .
docker push REGISTRY/pdvar:COMMIT
```

`workflow/containers/gauchian/Dockerfile` is the same pip pin (`gauchian==1.0.2`) if you want to skip the install at runtime. `workflow/containers/introme/DEPRECATED.md` records that Introme is the native WDL port, not a single user-built image.

## Internal sample JSON

The Assess task writes `sample.json` and passes it to `pdvar assess`. Keys: `sample_id`, `sex`, `vep_vcf`, `small_vcf`, `sv_vcf`, `sv_status`, `expansionhunter_vcf`, `gauchian_json`, `mtdna_vcf`, `introme_tsv`, `gene_bed`, `exon_bed`, `expression`, `splicevault`, `prs_weights`, `phenotype`, `qc`. This file is an in-task handoff. It is not an external manifest.

Introme's own tasks still accept `force_rerun_token` with default `""` so a cache-buster can be set on that subworkflow. `PdSample` does not expose it.
