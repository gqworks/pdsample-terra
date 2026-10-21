# PdSample workflow

WDL 1.0 for Terra / Cromwell. Entry workflow: `workflow/wdl/PdSample.wdl`.

Validate every WDL with:

```bash
bash workflow/scripts/validate_wdl.sh
```

The default path runs `womtool validate` on the committed files. PdSample imports:

- `workflow/wdl/vendor/gatk-sv/GATKSVPipelineSingleSample.wdl` (gatk-sv **v1.1.1** with absolute nested HTTP imports and HERVK/LINE1/`n_per_split` pass-throughs)
- `https://raw.githubusercontent.com/broadinstitute/gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl`

`miniwdl check` uses that same tree. It substitutes `stubs/MitochondriaPipeline.wdl` because upstream `AlignAndCall.wdl` multiplies an optional `Int?` and miniwdl rejects that. `gcp` is how the Introme GPU tasks pass `gpuCount`, `gpuType`, and `nvidiaDriverVersion` to Cromwell. miniwdl warns that it does not know the `gcp` key. That warning is expected and does not fail the check. The object is empty unless `use_gpu` is true, because Cromwell rejects `gpuCount: 0` while creating the job and never starts a VM.

`bash workflow/scripts/validate_wdl.sh --offline` replaces the vendored GATK-SV import and the mitochondria import with `workflow/wdl/stubs/`. It is not the Dockstore or Terra check.

## GATK-SV tag

`GATKSVPipelineSingleSample` is pinned to **v1.1.1** (Dockstore `SingleSamplePipeline:v1.1` / `v1.1.1`). Defaults match Broad: **Scramble on**, **MELT off**, new 1KG panel under `gatk-sv-ref-panel-1kg-v1-0-1`, clustering/GQ filtering inputs required.

Upstream v1.1.1 does not declare pass-throughs for nested `MakeCohortVcf.HERVK_reference`, `MakeCohortVcf.LINE1_reference`, and `RefineComplexVariants.n_per_split`. Cromwell therefore rejects an unmodified remote import as a sub-workflow. PdSample vendors a thin patched copy that wires those three inputs (defaults are the public HERVK/LINE1 beds and `n_per_split=15000` from Broad's test template). Nested WDLs still load from `raw.githubusercontent.com/.../v1.1.1/wdl/`.

GATK-SV is a required call. Panel and resource `File` inputs are required `PdSample` inputs; see `workflow/terra/gatksv_panel_keys.md` and `PdSample.example.json`. Leave `sv_*_docker` unset in Terra so the v1.1.1 pins apply.

## Task graph

| Step | Image | Version |
| --- | --- | --- |
| BAM header | `quay.io/biocontainers/samtools:1.24--h9dcdb79_1` | samtools 1.24 |
| Reference dict, HaplotypeCaller, hard filters | `us.gcr.io/broad-gatk/gatk:4.6.1.0` | GATK 4.6.1.0 |
| Reference check, QC summary | `python:3.11.11-slim-bookworm` | stdlib scripts |
| VerifyBamID2 | `quay.io/biocontainers/verifybamid2:2.0.3--hc004090_0` | 2.0.3 |
| somalier sex | `quay.io/biocontainers/somalier:0.3.5--h5205c93_0` | 0.3.5 |
| mosdepth | `quay.io/biocontainers/mosdepth:0.3.14--h87be163_2` | 0.3.14 |
| PASS and norm | `quay.io/biocontainers/bcftools:1.24--h118bc1c_4` | bcftools 1.24 |
| VEP | `quay.io/biocontainers/ensembl-vep:116.2--pl5321h2a3209d_0` | 116.2 |
| AlphaGenome AVI (precomputed table) | `python:3.11.11-slim-bookworm` | pip `pysam==0.24.1`, local Tabix lookup |
| ExpansionHunter | `quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5` | 5.0.0, 19-locus catalog |
| Gauchian | `python:3.11.11-slim-bookworm` unless `gauchian_docker` is set | pip `gauchian==1.0.2` |
| mtDNA | images inside MitochondriaPipeline | GATK 4.7.0.0 |
| GATK-SV | pins in `PdSample.wdl` (`sv_*_docker`) | gatk-sv v1.1.1, required |
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

GATK-SV image pins (overridable), from v1.1.1 `inputs/values/dockers.json`:

- manta `us.gcr.io/broad-dsde-methods/gatk-sv/manta:2023-09-14-v0.28.3-beta-3f22f94d`
- scramble `us.gcr.io/broad-dsde-methods/gatk-sv/scramble:2024-10-25-v0.29-beta-5ea22a52`
- wham / sv-base / sv-base-mini / samtools-cloud `2024-10-25-v0.29-beta-5ea22a52`
- sv-pipeline `2025-10-02-v1.1-483973d6`
- cnmops `2025-09-02-v1.0.5-f091af0b`
- gatk `us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT`
- gq recalibrator `us.gcr.io/broad-dsde-methods/markw/gatk:mw-tb-form-sv-filter-training-data-eb62a3`
- cloud SDK `google/cloud-sdk`

## Assess image

```bash
docker build -f workflow/containers/pdvar/Dockerfile -t REGISTRY/pdvar:COMMIT .
docker push REGISTRY/pdvar:COMMIT
```

Set `PdSample.pdvar_docker` to that image. Leave `pdvar_git_commit` empty when the image already contains `pdvar`.
