# Resources to stage

`PdSample` localizes these files. Example JSON paths are in `workflow/terra/PdSample.example.json`. Files that are not on public GCS, and the commands to copy them into a workspace bucket, are in [`workflow/terra/STAGING.md`](../workflow/terra/STAGING.md). `workflow/scripts/download_resources.sh` fetches the VerifyBamID files and the chr-prefixed somalier hg38 sites VCF (and the VEP cache when you pass `--vep`).

## You must supply

- **GCP requester-pays project.** Broad reference reads need `requester_pays_project` (`gsutil -u PROJECT`).
- **GRCh38 FASTA, index, and dictionary.** `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta` plus `.fai` and `.dict`. The CRAM must match this reference. The header check fails on an MD5 or `chr` prefix mismatch.
- **VEP 116 cache, as a `.tar.gz`.** `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz` (about 26 GB). The VEP task extracts it. There is no public GCS copy. Stage it with the commands in [`workflow/terra/STAGING.md`](../workflow/terra/STAGING.md). The example path is `gs://YOUR_WORKSPACE_BUCKET/resources/vep/homo_sapiens_vep_116_GRCh38.tar.gz`.
- **VerifyBamID2 2.0.3 SVD files (four inputs).** Public raw URLs, also fetched by `download_resources.sh`:
  - `https://github.com/Griffan/VerifyBamID/raw/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.UD`
  - `https://github.com/Griffan/VerifyBamID/raw/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.V`
  - `https://github.com/Griffan/VerifyBamID/raw/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.bed`
  - `https://github.com/Griffan/VerifyBamID/raw/v2.0.3/resource/1000g.phase3.100k.b38.vcf.gz.dat.mu`
- **somalier sites VCF.** Chr-prefixed hg38 sites linked from the somalier v0.3.5 release notes: `https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz` (about 260 KB). The v0.3.5 release asset list is only the binary. Do not use `sites.hg38.nochr.vcf.gz`. Stage it with the commands in [`workflow/terra/STAGING.md`](../workflow/terra/STAGING.md). The example path is `gs://YOUR_WORKSPACE_BUCKET/resources/somalier/sites.hg38.vcf.gz`.
- **mtDNA bundle for MitochondriaPipeline 4.7.0.0.** Required inputs are `mt_fasta`, `mt_fasta_index`, `mt_dict`, `mt_amb`, `mt_ann`, `mt_bwt`, `mt_pac`, `mt_sa`, `blacklisted_sites`, `blacklisted_sites_index`, the shifted FASTA/index/dict/BWA indexes, `shift_back_chain`, `control_region_shifted_reference_interval_list`, and `non_control_region_interval_list`. Both JSON files point at the public hg38 chrM bundle under `gs://gcp-public-data--broad-references/hg38/v0/chrM/`, the same objects as GATK 4.7.0.0 `ExampleInputsMitochondriaPipeline.json`. Contig name is `chrM`. `blacklisted_sites` is `blacklist_sites.hg38.chrM.bed` with `blacklist_sites.hg38.chrM.bed.idx` beside it. That tag's `FilterMutectCalls` `--mask` takes this BED. The public prefix has no blacklist VCF. The same `requester_pays_project` used for the genome FASTA covers this bucket.
- **Introme hg38 GTF.** Example uses GENCODE v44: `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz`. The Introme pin is `d6da48bc84abda66f588470f79abcf1244acfb21`.
- **`pdvar_git_commit` or a custom image.** The example uses `python:3.11.11-slim-bookworm` and `pdvar_git_commit` `main`, which installs `git+https://github.com/gqworks/pdsample-terra.git@main`. After the public snapshot exists, pin `pdvar_git_commit` to that commit SHA of `gqworks/pdsample-terra`. Or build `workflow/containers/pdvar/Dockerfile` and set `pdvar_docker`.
- **GATK-SV v0.29-beta reference panel and hg38 resources.** Structural-variant calling is required. `workflow/terra/PdSample.example.json` already points at the public files tag v0.29-beta resolves from `inputs/values/ref_panel_1kg.json` and `inputs/values/resources_hg38.json` (`build_inputs.py` with the single-sample test template and `ref_panel_1kg`). `workflow/terra/PdSample.terra_inputs.json` binds the same inputs to workspace attributes. The full table is `workflow/terra/gatksv_panel_keys.md`. v1.1.1 is not used: its single-sample workflow cannot be called as a sub-workflow.

### One matched panel

Every panel file has to come from one bundle built with a GATK-SV version compatible with v0.29-beta. Do not mix a cohort `ref_panel_vcf` or ped with the 1KG gCNV models or bincov matrix. The files that must move together are `ref_panel_vcf`, `ref_panel_ped_file`, `ref_panel_samples_list`, `ref_panel_bincov_matrix`, the PE/SR/SD file lists, `ref_panel_contig_ploidy_model_tar`, `ref_panel_gcnv_model_tars_list`, the std Manta and Wham VCF tars, the del and dup beds, `ref_panel_cutoffs`, the four genotype separation cutoffs, and the PE/SR metrics.

The public `ref_panel_1kg` JSON already combines more than one path family, and that combination is what the tag ships:

- `contig_ploidy_model_tar` is the 1KG **v2** gCNV tar under `gs://gatk-sv-resources-public/.../1KG/v2/gcnv/`. The gCNV model shards named by that bundle are the same public v2 tars. The list file, sample list, bincov matrix, del/dup beds, cutoffs, genotype separation cutoffs, and PE/SR metrics are under batch `38c65ca4-2a07-4805-86b6-214696075fef`.
- std Manta and Wham tars are GatherBatchEvidence `b705e9f6-c401-4a9a-a059-1a68d1ca5f98`.
- PE/SR/SD file lists are `tws-no-cram-conversion/GatherSampleEvidenceBatch`.
- The cleaned VCF is `mw-vcf-reshard/ref_panel_1kg.cleaned.vcf.gz`. The ped is public v1 `1kg_ref_panel_v1.ped`.
- `ref_panel_1kg.json` also has a batch `qc_definitions` TSV. The single-sample input is not that file. The Terra template sets `workspace.ref_panel_qc_definitions` from `resources_hg38.single_sample_qc_definitions` (`.../1KG/v2/single_sample.qc_definitions.tsv`).

Keep that set intact. Replacing only the ploidy tar, or only the VCF, breaks the match.

### Your own Terra panel

If the workspace was created from the GATK-SV featured workspace, the attributes are already there. Import `workflow/terra/PdSample.terra_inputs.json`, or set each input in the Terra Inputs tab to the expression in `gatksv_panel_keys.md`. Sample columns stay `this.sample_id`, `this.cram`, `this.crai`, and `this.sex`.

Panel attributes use the `ref_panel_*` names (`workspace.ref_panel_vcf`, `workspace.ref_panel_ped_file`, `workspace.ref_panel_samples_list`, and the rest of that family). hg38 resource files use the `reference_*` twins (`workspace.reference_fasta`, `workspace.reference_genome_file`, `workspace.reference_empty_file`, and so on). `pesr_exclude_intervals` reads `workspace.reference_pesr_exclude_list`. `qc_definitions` reads `workspace.ref_panel_qc_definitions`.

Docker inputs are not workspace attributes in this method config. Featured-workspace docker values can mix v0.29-beta images with other tags (`sv_pipeline_docker` / `sv_pipeline_qc_docker` in particular). Leave `sv_*_docker` unset so `PdSample.wdl` keeps the v0.29-beta pins.

## Optional

- **gnomAD and ClinVar VCFs** for VEP `--custom`, plus their tabix indexes (`gnomad_vcf_index`, `clinvar_vcf_index`). Without gnomAD, the Introme rare filter has no allele frequency to apply. Missing AF is kept by the assessor.
- **VEP plugin data** (CADD, REVEL, SpliceAI) via `vep_extra_args` if you have a license and the files. The scorer reads those CSQ fields when they are present.
- **Introme assets tar.** Optional `introme_assets_tar`: a workspace copy of `https://github.com/CCICB/introme/archive/d6da48bc84abda66f588470f79abcf1244acfb21.tar.gz`. Unset, the fetch task downloads that pin.
- **Introme image mirrors.** Defaults are the Docker Hub tags in `workflow/wdl/Introme.wdl`. Mirror them to Artifact Registry and override the `*_docker` inputs if Docker Hub pulls are blocked. SpliceAI, Pangolin, and MMSplice have CPU and GPU tags. `introme_use_gpu` selects the GPU tags.
- **Gauchian image.** `workflow/containers/gauchian/Dockerfile` installs `gauchian==1.0.2`. The default task pip-installs that pin into `python:3.11.11-slim-bookworm`.
- **Expression, SpliceVault, PRS weights.** File inputs on Assess. Expression and PRS numbers in config stay null until cited, so supplying a table does not by itself change the score.
- **`age_onset` and `symptoms`.** Optional strings copied into the sample JSON phenotype block.

## Already in the repo

- Calling BED: `workflow/assets/monopdaus-genes.grch38.bed`
- SV overlap BEDs: `resources/panel/monopdaus-genes.grch38.bed`, `resources/panel/monopdaus-exons.grch38.bed`
- Inheritance table: `resources/panel/monopdaus-gene-inheritance.csv` (placeholder labels)
- ExpansionHunter catalog: `workflow/assets/expansionhunter_pd.hg38.json` (19 loci; NOTCH2NLC, FGF14, BEAN1, and LRP12 are not in this catalog)
- Config: `config/config.yaml`, `config/columns.yaml`, `config/gene_phenotypes.yaml`
- QC scripts: `workflow/scripts/check_reference.py`, `workflow/scripts/qc_summary.py`

Panel inheritance labels, penetrance, PRS beta, expression cutoffs, and STR pathogenic repeat counts are null or placeholder until you cite them. The workflow will still rank variants.
