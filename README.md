# Per-sample Parkinson's variant assessment

One Terra sample row goes in as a BAM or CRAM. `PdSample` calls every required tool from that alignment, annotates the small-variant VCF, and writes a ranked candidate list for that person. A cohort is a Terra `sample_set`: Terra launches `PdSample` once per row. There is no cohort workflow.

```mermaid
flowchart TD
  cram[BAM or CRAM plus reference] --> check[Reference check]
  check --> qc[VerifyBamID2, somalier sex, mosdepth]
  qc --> hc[HaplotypeCaller GVCF]
  hc --> gt[GenotypeGVCFs and hard filters]
  gt --> norm[PASS and bcftools norm]
  norm --> vep[VEP 116]
  vep --> lift[Lift gnomAD INFO]
  lift --> introme[Introme v2 on the panel BED]
  norm --> avi[AlphaGenome Atlas AVI]
  qc --> eh[ExpansionHunter]
  qc --> gba[Gauchian]
  qc --> mito[MitochondriaPipeline]
  qc --> gatksv[GATKSVPipelineSingleSample]
  pheno[Phenotype report or age, symptoms, ancestry] --> assess[Assess]
  introme --> assess
  avi --> assess
  eh --> assess
  gba --> assess
  mito --> assess
  gatksv --> assess
  assess --> row[Candidates, report, QC, callsets]
```

## Dockstore

Register the repository root. `.dockstore.yml` publishes one workflow:

| Name | Descriptor | Example JSON |
| --- | --- | --- |
| `PdSample` | `/workflow/wdl/PdSample.wdl` | `/workflow/terra/PdSample.example.json` |

Import `PdSample` into a Terra workspace. Bind `sample_id`, `bam_or_cram_file`, `bai`, `sex`, `age_onset`, `symptoms`, `ethnicity`, and `phenotype_report` to the sample table (method inputs are still named `cram` and `crai` in the WDL). Put reference and resource files on the workspace data table or in the method configuration.

Broad hg38 paths under `gs://gcp-public-data--broad-references` are requester pays. On submit, enable **Use reference disks** and set workflow option `google_project` to your billing GCP project (`workflow/terra/PdSample.workflow_options.example.json`). Set `PdSample.requester_pays_project` to that same project. See [`docs/setup_checklist.md`](docs/setup_checklist.md).

Terra GCP Batch only localizes `gs://` `File` inputs. Relative paths and `https://` URLs fail if they are workflow `File` values, including WDL defaults. `FetchRepoAssets` downloads panel, config, scripts, and VerifyBamID SVD files from GitHub inside the task (`pdvar_git_commit`, default `main`) and passes the resulting `gs://` outputs to later tasks. Stage VEP, somalier, GENCODE GTF, and the alignment into your workspace bucket per [`workflow/terra/STAGING.md`](workflow/terra/STAGING.md). Re-import `PdSample.terra_inputs.json` after upgrades and delete leftover `../` or `https://` method inputs.

## Terra sample table

| Column | Required | Meaning |
| --- | --- | --- |
| `sample_id` | yes | Identifier written into outputs. Use a synthetic id in examples. |
| `bam_or_cram_file` | yes | Aligned GRCh38 BAM or CRAM for this row. |
| `bai` | yes | BAM index (`.bai`) or CRAM index (`.crai`). |
| `sex` | no | `male` or `female`. Blank, `unspecified`, or anything else is inferred from somalier, then from the X heterozygosity rate. |
| `age_onset` | no | Age at onset in years. Blank skips onset concordance. A numeric string is accepted. |
| `symptoms` | no | Free-text clinical features. Keywords are matched to each gene profile in `config/gene_phenotypes.yaml` (for example dementia, tremor, dystonia, levodopa responsive, hallucinations, autonomic). |
| `ethnicity` | no | Ancestry text. Matched against `enriched_ancestry` on the gene profile (for example Ashkenazi Jewish for GBA1 and LRRK2). |
| `phenotype_report` | no | `gs://` text file in the GP2-style sheet layout (`GP2ID`, sex, age at onset, levodopa, exposures and comorbidities). Leave the cell empty when the patient has no report. |

Add `age_onset`, `symptoms`, `ethnicity`, and `phenotype_report` to the sample table even when a row leaves them blank. `workflow/terra/PdSample.terra_inputs.json` binds those columns. A phenotype report is parsed into age at onset, levodopa, sampling age, diagnosis, and the comorbidity list. Explicit `age_onset`, `symptoms`, and `ethnicity` override the same fields from the report. The same text pasted into `symptoms` is recognized and parsed instead of being searched as free-text keywords.

Assess writes `clinical_match` on each candidate, including risk variants (`risk_carrier`, `risk_biallelic`, and genes flagged `risk_allele_exception`). Onset, symptom keywords, and ancestry drive that match. Levodopa use is stored as a biomarker note on early-onset recessive genes. Comorbidities are printed on the report. They are not gene-profile keywords, so hypertension, constipation, hearing loss, and the other lines on the example sheet do not by themselves mark a risk variant as clinically matched. Genes with no profile stay `no_profile`. A blank phenotype leaves `clinical_match` as `not_assessed`.

A `sample_set` is only a list of those rows. Terra runs `PdSample` per row. Shared reference inputs stay in the method config.

`workflow/terra/PdSample.example.json` uses `SYNTHETIC_SAMPLE`, `gs://YOUR_WORKSPACE_BUCKET/...` placeholders for the example BAM and for the VEP cache and somalier sites (see [`workflow/terra/STAGING.md`](workflow/terra/STAGING.md)), raw GitHub URLs on `main` for files in this repo, the public VerifyBamID 2.0.3 resource files, the public GATK hg38 chrM bundle, and the public GATK-SV v0.29-beta hg38 resources plus the 1KG reference panel. `pdvar_git_commit` defaults to `main` on https://github.com/gqworks/pdsample-terra. After this snapshot is published, pin that input to the commit SHA of `pdsample-terra` so Terra does not float. `workflow/terra/PdSample.terra_inputs.json` is the Terra form of that config: sample fields are `this.sample_id`, `this.bam_or_cram_file`, `this.bai`, `this.sex`, `this.age_onset`, `this.symptoms`, `this.ethnicity`, and `this.phenotype_report`, the mitochondria files are the same public `gs://` paths, `vep_cache` and `somalier_sites` stay `gs://YOUR_WORKSPACE_BUCKET/...`, and each GATK-SV file is a `workspace.<key>` expression. The mapping is `workflow/terra/gatksv_panel_keys.md`.

## What each stage does

1. **Reference check.** The alignment header `@SQ` MD5s are compared with the reference dictionary. A missing `chr` prefix fails the check.
2. **QC.** VerifyBamID2 contamination (`FREEMIX` above 0.02 fails the sample by default), somalier extract plus relate on this sample only (sex; relatedness is not run), and mosdepth on the calling BED. Callers wait on the QC gate.
3. **Small variants.** GATK 4.6.1.0 HaplotypeCaller in GVCF mode, then GenotypeGVCFs, then GATK hard filters, then `bcftools view -f PASS` and `bcftools norm -m -both`. Intervals default to `workflow/assets/monopdaus-genes.grch38.bed`. Set `wgs_small_variants` true to call the whole genome.
4. **Structural variants.** `GATKSVPipelineSingleSample` v0.29-beta runs on every sample. The reference-panel and hg38 resource inputs are required. `sv_status` is `called`. A missing SV VCF is an error in `pdvar assess`.
5. **STRs.** ExpansionHunter 5.0.0 with the vendored 19-locus catalog `workflow/assets/expansionhunter_pd.hg38.json`.
6. **GBA1.** Gauchian 1.0.2.
7. **mtDNA.** GATK `MitochondriaPipeline` 4.7.0.0.
8. **VEP 116.2.** Offline cache, GRCh38, `--pick`. Optional `--custom` gnomAD (`AF|AF_popmax|nhomalt`) and ClinVar. A following task copies `INFO/gnomAD` into `INFO/gnomAD_AF` and `INFO/AF_popmax` so the Introme rare filter can see them.
9. **AlphaGenome Atlas AVI.** PASS SNVs are looked up in the precomputed AlphaGenome AVI SNV table (Tabix, GRCh38) when `avi_table` is set. No API key or network call is needed. Phred 30 (top 0.1% of genome-wide SNVs) and Phred 20 (top 1%) add triage points and break rank ties. An unset `avi_table` leaves the module `not_assessed`. Indels are not scored. This is a research rank, not a diagnostic call.
10. **Introme v2.** Native WDL port of CCICB/introme `d6da48bc84abda66f588470f79abcf1244acfb21` (`v2@d6da48b`). Restricted to the panel BED. SpliceAI, MMSplice, Pangolin, SPiP, and Spliceogen run in parallel after VariantInfo.
11. **Assess.** `pdvar assess` ranks candidates and writes the Terra outputs. A phenotype report, or age at onset, symptoms, and ancestry, is compared with `config/gene_phenotypes.yaml`. `clinical_match` on a risk variant records whether that gene's onset range, key clinical features, or enriched ancestry fit the patient. The match is a triage note. It does not add evidence points. A higher `phenotype_score` breaks ties after tier and AlphaGenome Phred. Profiles cover every gene in `resources/panel/monopdaus-gene-inheritance.csv` (sources and confidence in `config/gene_phenotypes.SOURCES.md`). They are a starter literature synthesis and need expert clinical review. Low-evidence candidates and the risk genes PSAP and TMEM175 have minimal profiles with no onset range, so they are usually `uninformative`. A gene outside the panel is `no_profile`.

## Hard filters

Single-sample hard filters are used because VQSR fits tranches on a cohort and this repository does not ship a CNN scoring model. SNP sites fail when `QD < 2.0`, `QUAL < 30.0`, `SOR > 3.0`, `FS > 60.0`, `MQ < 40.0`, `MQRankSum < -12.5`, or `ReadPosRankSum < -8.0`. INDEL sites fail when `QD < 2.0`, `QUAL < 30.0`, `FS > 200.0`, or `ReadPosRankSum < -20.0`. Only `PASS` records continue.

## Outputs written back to the row

| Output | File |
| --- | --- |
| `candidates_tsv` | Ranked variants that are not `not_qualifying` |
| `report_md` / `report_html` | Per-sample triage report |
| `qc_summary` | QC JSON, including module status and `sv_status` |
| `liability_tsv` / `followup_tsv` | Explained-liability row and single-hit follow-up |
| `gvcf`, `small_vcf` | HaplotypeCaller GVCF and filtered small-variant VCF |
| `vep_vcf` | VEP VCF after the gnomAD INFO lift |
| `sv_vcf` / `sv_vcf_index` | GATK-SV single-sample VCF and index |
| `sv_status` | `called` |
| `expansionhunter_vcf`, `gauchian_json`, `mtdna_vcf` | STR, GBA1, and mtDNA callsets |
| `introme_tsv`, `introme_vcf` | Introme scores |
| `sex_used` | Provided sex, or the inferred value |
| `sample_json` | Internal JSON the Assess task reads |

## Assess container

Build from the repository root and push:

```bash
docker build -f workflow/containers/pdvar/Dockerfile -t REGISTRY/pdvar:COMMIT .
docker push REGISTRY/pdvar:COMMIT
```

Set `pdvar_docker` to that image. The image copies this package, `config/`, and `resources/panel/`.

To run without a custom image, leave `pdvar_docker` at `python:3.11.11-slim-bookworm`. `pdvar_git_commit` defaults to `main`, and the task installs `git+https://github.com/gqworks/pdsample-terra.git@main`. Pin `pdvar_git_commit` to a commit SHA of `gqworks/pdsample-terra` for a reproducible run. An empty commit with no installed `pdvar` stops the task.

## Two panel BEDs

`workflow/assets/monopdaus-genes.grch38.bed` is the calling and Introme interval list (gene spans). `resources/panel/monopdaus-genes.grch38.bed` and `monopdaus-exons.grch38.bed` are the overlap BEDs the SV scorer uses. Inheritance labels live in `resources/panel/monopdaus-gene-inheritance.csv` and are placeholders pending expert review. `GBA` is aliased to `GBA1`.

## Cost

The default interval list keeps HaplotypeCaller, VEP, and Introme on the panel gene spans. Whole-genome small-variant calling (`wgs_small_variants`) and a full GATK-SV reference-panel run dominate cost. Introme's SpliceAI, MMSplice, and Pangolin tasks are the next largest CPU (or GPU) block. ExpansionHunter, Gauchian, VerifyBamID2, somalier, and mosdepth are small beside those. The mitochondria workflow realigns chrM. There is no joint genotyping and no VQSR.

`ReleaseAfterQc` hardlinks the CRAM after the QC gate so the imported mitochondria and GATK-SV workflows do not start early. The task disk is twice the CRAM size so a cross-filesystem copy can finish when a hardlink is unavailable.

## Local assessment

```bash
pip install -e ".[dev]"
pdvar assess sample.json --panel resources/panel/monopdaus-gene-inheritance.csv --outdir out
pytest
bash workflow/scripts/validate_wdl.sh
```

`pdvar expression-urls` prints GTEx and gnomAD pext locations and downloads nothing. `scripts/fetch_expression_resources.py` builds an optional expression table.

## Scientific scope

The scorer keeps panel genes, ACMG-style point evidence (ClinVar P/LP, predicted null alleles, CADD, REVEL, SpliceAI, absent gnomAD, Introme), a second allele inside the same sample, SV/STR/GBA1/mtDNA rows, isoform and expression weights, phenotype concordance, and penetrance scaffolding. Point totals are a triage rank. They are not an ACMG/AMP classification. Phenotype concordance asks whether a candidate gene, including a risk variant, fits the patient's age at onset, symptoms, and ancestry. Penetrance, PRS coefficients, expression cutoffs, and STR pathogenic repeat counts stay null until a citation is set in `config/config.yaml`. TTN truncations are tagged `incidental` and are not relabelled benign. Compound hets stay `AR_comphet_candidate` unless the VCF already carries a shared `PS` tag.

Resources you still have to stage are listed in [`workflow/terra/STAGING.md`](workflow/terra/STAGING.md) and [`docs/setup_checklist.md`](docs/setup_checklist.md). Tool versions and images are in [`workflow/README.md`](workflow/README.md).
