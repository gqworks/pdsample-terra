# Files to stage into the workspace bucket

`PdSample.example.json` and `PdSample.terra_inputs.json` point at public `gs://` objects wherever an anonymous copy exists. The inputs below do not. Their values stay `gs://YOUR_WORKSPACE_BUCKET/...`. Replace `YOUR_WORKSPACE_BUCKET` with the Terra workspace bucket when you launch. Do not commit a real bucket name, sample id, or billing project.

Mitochondria references, the GRCh38 FASTA, VerifyBamID 2.0.3 SVD files (vendored under `workflow/assets/verifybamid/`), panel/config/scripts bundled with the Dockstore workflow, and the GATK-SV v1.1.1 public resources are not staged here. Terra GCP Batch only localizes `gs://` (and DRS) `File` inputs, not `https://` URLs.

See [`docs/setup_checklist.md`](../../docs/setup_checklist.md).

## Placeholders

| Input | JSON | Placeholder | Upstream | Approx. size |
| --- | --- | --- | --- | --- |
| `PdSample.bam` | example only | `gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.bam` | your own GRCh38 BAM | the alignment |
| `PdSample.bai` | example only | `gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.bam.bai` | `.bai` for that BAM | small |
| `PdSample.vep_cache` | both | `gs://YOUR_WORKSPACE_BUCKET/resources/vep/homo_sapiens_vep_116_GRCh38.tar.gz` | Ensembl VEP cache release 116 | 26 GB (27,644,657,162 bytes) |
| `PdSample.somalier_sites` | both | `gs://YOUR_WORKSPACE_BUCKET/resources/somalier/sites.hg38.vcf.gz` | somalier v0.3.5 release notes, chr-prefixed hg38 | 260 KB (265,818 bytes) |
| `PdSample.introme_gtf` | both | `gs://YOUR_WORKSPACE_BUCKET/resources/gencode/gencode.v44.annotation.gtf.gz` | GENCODE v44 annotation | ~47 MB |
| `PdSample.avi_table` | terra only | `gs://YOUR_WORKSPACE_BUCKET/resources/alphagenome/alphagenome_variant_impact_score_snvs.tsv.gz` | AlphaGenome "AVI SNV scores" download, GRCh38, Tabix | 88.5 GB |
| `PdSample.avi_table_index` | terra only | `gs://YOUR_WORKSPACE_BUCKET/resources/alphagenome/alphagenome_variant_impact_score_snvs.tsv.gz.tbi` | Tabix index for that table | small |

`PdSample.terra_inputs.json` binds the alignment to `${this.bam_or_cram_file}` and `${this.bai}` on the sample table (WDL inputs are `bam` and `bai`). Only BAM + `.bai` is supported. The example JSON paths exist so Dockstore and `womtool` have a `File`. `SYNTHETIC_SAMPLE` is not a real sample.

`PdSample.requester_pays_project` stays `YOUR_GCP_PROJECT`. It is the GCP project that pays for requester-pays reads of `gs://gcp-public-data--broad-references` (the genome FASTA and the chrM bundle). It is not a file.

## Sample BAM

Copy a BAM aligned to `Homo_sapiens_assembly38` (chr-prefixed, dictionary MD5s matching that FASTA) and its `.bai` into the workspace bucket, then set the sample-table columns `bam_or_cram_file` and `bai`. CRAM is not supported. This repository does not publish an example alignment.

```bash
gcloud storage cp /path/to/SYNTHETIC_SAMPLE.bam \
  gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.bam
gcloud storage cp /path/to/SYNTHETIC_SAMPLE.bam.bai \
  gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.bam.bai
```

## VEP 116 cache

The VEP task runs `quay.io/biocontainers/ensembl-vep:116.2--pl5321h2a3209d_0` and extracts this tarball as `--dir_cache`. The cache release matches the Ensembl release in that image (116), not a later cache.

Source: `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz`

Checked 2026-10-09: HTTP 200, `Content-Length` 27,644,657,162, `Last-Modified` Thu, 09 Apr 2026 00:45:28 GMT. No copy was found under `gs://gcp-public-data--broad-references/hg38/v0/vep/`.

The download needs about 26 GB of free disk.

```bash
curl -fL -o homo_sapiens_vep_116_GRCh38.tar.gz \
  "https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz"
gcloud storage cp homo_sapiens_vep_116_GRCh38.tar.gz \
  gs://YOUR_WORKSPACE_BUCKET/resources/vep/homo_sapiens_vep_116_GRCh38.tar.gz
```

`workflow/scripts/download_resources.sh DEST --vep` fetches the same tarball into `DEST/vep/`. Copy that file with the `gcloud storage cp` line above.

## GENCODE v44 GTF (Introme)

Terra cannot read `https://ftp.ebi.ac.uk/...` as a workflow `File`. Copy GENCODE v44 into the workspace bucket once:

Source: `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz`

```bash
curl -fL -o gencode.v44.annotation.gtf.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
gcloud storage cp gencode.v44.annotation.gtf.gz \
  gs://YOUR_WORKSPACE_BUCKET/resources/gencode/gencode.v44.annotation.gtf.gz
```

Set `PdSample.introme_gtf` to that `gs://` path (already in `PdSample.terra_inputs.json`).

## somalier sites (hg38, chr-prefixed)

The somalier task runs `quay.io/biocontainers/somalier:0.3.5--h5205c93_0`. Release [v0.3.5](https://github.com/brentp/somalier/releases/tag/v0.3.5) ships the binary only. The release notes link the sites VCFs. Use the chr-prefixed hg38 file. The reference in this workflow is chr-prefixed, so `sites.hg38.nochr.vcf.gz` is the wrong file.

Source: `https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz`

Checked 2026-10-09: HTTP 200 after the GitHub redirect, 265,818 bytes, gzip, contigs `chr1`–`chr22`, `chrX`, `chrY`. No `.tbi` is published next to that file, and `somalier_sites` is a single `File` (the task does not take an index). No copy was found under `gs://gcp-public-data--broad-references/hg38/v0/` with a `somalier` prefix.

```bash
curl -fL -o sites.hg38.vcf.gz \
  "https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz"
gcloud storage cp sites.hg38.vcf.gz \
  gs://YOUR_WORKSPACE_BUCKET/resources/somalier/sites.hg38.vcf.gz
```

`workflow/scripts/download_resources.sh` also downloads this file into `DEST/somalier/`.

## AlphaGenome AVI SNV table

`AnnotateAvi` looks up PASS SNVs in the precomputed AlphaGenome AVI table with a local Tabix query (`pysam==0.24.1`). It does not call the AlphaGenome Atlas API, so it needs no key and works in workspaces where that API returns HTTP 403 (for example AnVIL).

Source: https://alphagenome.google/downloads, "AVI SNV scores" (Tabix, 88.5 GB, permissive use). The download needs a signed-in browser, so it cannot be fetched with `curl` here. Download it once, then copy the table and its index into the bucket:

```bash
gcloud storage cp alphagenome_variant_impact_score_snvs.tsv.gz \
  gs://YOUR_WORKSPACE_BUCKET/resources/alphagenome/alphagenome_variant_impact_score_snvs.tsv.gz
gcloud storage cp alphagenome_variant_impact_score_snvs.tsv.gz.tbi \
  gs://YOUR_WORKSPACE_BUCKET/resources/alphagenome/alphagenome_variant_impact_score_snvs.tsv.gz.tbi
```

Use the file names you actually downloaded and set `PdSample.avi_table` and `PdSample.avi_table_index` to match. If no index came with the table, build one with `tabix -s 1 -b 2 -e 2 alphagenome_variant_impact_score_snvs.tsv.gz` (the file must be bgzip, not plain gzip). Check the header first:

```bash
gzip -dc alphagenome_variant_impact_score_snvs.tsv.gz | head -3
```

The script expects `CHROM POS REF ALT raw_score PHRED` and also accepts common alternatives (for example `avi_score`, `phred_score`). `chr1` and `1` contig names both work. Cromwell copies the whole table into each `AnnotateAvi` task, and the task disk is sized from the table. Leave both inputs unset to skip AVI (`avi: not_assessed`).

## Terra notebook

Terra notebooks set `WORKSPACE_BUCKET` to this workspace's bucket. The VEP cache is about 26 GB.

```python
import os
import subprocess

bucket = os.environ["WORKSPACE_BUCKET"]

def stage(url: str, object_name: str) -> None:
    local = "/tmp/" + object_name.rsplit("/", 1)[-1]
    subprocess.check_call(["curl", "-fL", "-o", local, url])
    subprocess.check_call(["gcloud", "storage", "cp", local, f"{bucket}/{object_name}"])

stage(
    "https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz",
    "resources/vep/homo_sapiens_vep_116_GRCh38.tar.gz",
)
stage(
    "https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz",
    "resources/somalier/sites.hg38.vcf.gz",
)
```

Point `PdSample.vep_cache` and `PdSample.somalier_sites` at those two objects. Put the alignment and index on the sample table (`bam_or_cram_file`, `bai`).

## VEP pathogenicity annotations (academic dbNSFP 5.4a)

gnomAD AF, ClinVar significance, CADD, REVEL, and the other in-silico scores all come from **dbNSFP 5.4a** (no separate gnomAD / ClinVar VCF `--custom` inputs). SpliceAI for triage is from **Introme**, not VEP.

| Input | Staged path | Source |
| --- | --- | --- |
| `dbnsfp_variant_chroms` | `gs://YOUR_WORKSPACE_BUCKET/resources/dbnsfp/dbNSFP5.4a_variant.chr{1-22,X,Y}.gz` | Academic **dbNSFP 5.4a** as distributed (plain gzip per chromosome). `PrepareDbnsfpForVep` builds a panel-subset VEP-ready bgzip+tabix. Optional `chrM` can be added if needed. |
| `dbnsfp` / `dbnsfp_index` (optional) | `gs://…/dbNSFP5.4a_grch38.gz` (+ `.tbi`) | Pre-built VEP-ready file; when set, skips per-chrom prep |

### Stage academic dbNSFP 5.4a (per-chromosome gz)

1. Download the academic (`a`) branch from [dbNSFP releases](https://dbnsfp.org) (registration required; CC BY-NC-ND + upstream non-commercial terms for CADD/REVEL and others).
2. Upload the distributed files as-is: `dbNSFP5.4a_variant.chr1.gz` … `chr22.gz`, `chrX.gz`, `chrY.gz` (and `chrM.gz` if desired) under `resources/dbnsfp/`. Do **not** rename or concatenate offline unless you prefer a pre-built VEP-ready file.
3. Point `PdSample.dbnsfp_variant_chroms` at that array (pre-filled in `PdSample.terra_inputs.json`). The workflow’s `PrepareDbnsfpForVep` task subsets to the panel BED, then bgzips and tabixes for the VEP plugin.
4. Default `dbnsfp_fields` matches `dbNSFP5.4a_variant.columns.tsv`: SIFT/PolyPhen2, Eigen, MutationTaster, fathmm-XF, `CADD_phred`, `REVEL_score`, `clinvar_clnsig`, `gnomAD4.1_joint_AF`, `gnomAD4.1_joint_POPMAX_AF`. Classic `FATHMM_*` is **not** in 5.4a.

SIFT and PolyPhen also come from the VEP cache via `--sift b --polyphen b` (`vep_sift_polyphen`, default true). `LiftGnomadInfo` copies `gnomAD4.1_joint_*` from CSQ into `INFO/gnomAD_AF` / `INFO/AF_popmax` for Introme’s rare filter. `pdvar` maps `clinvar_clnsig` → `clinvar_class`.

## Optional inputs that are not in the JSON

These are unset in both JSON files. A default launch does not need them.

- **Introme GTF.** Stage GENCODE v44 to `gs://YOUR_WORKSPACE_BUCKET/resources/gencode/gencode.v44.annotation.gtf.gz` (see above). The Introme pin is `d6da48bc84abda66f588470f79abcf1244acfb21`.
- **Introme assets tar.** Optional `introme_assets_tar`. Unset, the fetch task downloads `https://github.com/CCICB/introme/archive/d6da48bc84abda66f588470f79abcf1244acfb21.tar.gz`.
- **Genome-wide CADD / REVEL plugins** (`cadd_snv`, optional `cadd_indels`, `revel`, and indexes). Optional fallback when you need scores outside dbNSFP’s missense/splice SNV scope; academic 5.4a already provides `CADD_phred` / `REVEL_score` for those sites.
- **ExpansionHunter / Gauchian catalogs.** `eh_catalog` is already the vendored 19-locus hg38 catalog in this repo (ExpansionHunter 5.0.0). Gauchian 1.0.2 has no catalog file input.
