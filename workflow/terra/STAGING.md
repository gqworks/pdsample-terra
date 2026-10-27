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

## VEP frequency and pathogenicity annotations

`VepAnnotate` wires these as first-class `File` inputs (not only `vep_extra_args`). Stage each file **and its tabix index** under the workspace bucket, then set the matching `PdSample.*` input (see `PdSample.terra_inputs.json`). The example JSON leaves them unset so a synthetic dry-run does not localize multi-GB resources.

| Input | Staged path | Source |
| --- | --- | --- |
| `gnomad_genomes_vcf` | `gs://YOUR_WORKSPACE_BUCKET/resources/gnomad/gnomad.genomes.v4.1.sites.vcf.gz` | gnomAD v4.1 genomes sites (GRCh38); subset to the panel BED if full sites are too large |
| `gnomad_genomes_vcf_index` | `gs://YOUR_WORKSPACE_BUCKET/resources/gnomad/gnomad.genomes.v4.1.sites.vcf.gz.tbi` | tabix index for genomes sites |
| `gnomad_exomes_vcf` | `gs://YOUR_WORKSPACE_BUCKET/resources/gnomad/gnomad.exomes.v4.1.sites.vcf.gz` | gnomAD v4.1 exomes sites (GRCh38) |
| `gnomad_exomes_vcf_index` | `gs://YOUR_WORKSPACE_BUCKET/resources/gnomad/gnomad.exomes.v4.1.sites.vcf.gz.tbi` | tabix index for exomes sites |
| `clinvar_vcf` | `gs://YOUR_WORKSPACE_BUCKET/resources/clinvar/clinvar.vcf.gz` | NCBI ClinVar GRCh38 VCF |
| `clinvar_vcf_index` | `gs://YOUR_WORKSPACE_BUCKET/resources/clinvar/clinvar.vcf.gz.tbi` | tabix index for ClinVar |
| `cadd_snv` | `gs://YOUR_WORKSPACE_BUCKET/resources/cadd/whole_genome_SNVs.tsv.gz` | CADD (license required); VEP 116 `--plugin CADD,snv=...` |
| `cadd_snv_index` | `gs://YOUR_WORKSPACE_BUCKET/resources/cadd/whole_genome_SNVs.tsv.gz.tbi` | tabix index for CADD SNVs |
| `cadd_indels` | `gs://YOUR_WORKSPACE_BUCKET/resources/cadd/gnomad.genomes.r4.0.indel.tsv.gz` | CADD indel scores (optional but recommended with `cadd_snv`) |
| `cadd_indels_index` | `gs://YOUR_WORKSPACE_BUCKET/resources/cadd/gnomad.genomes.r4.0.indel.tsv.gz.tbi` | tabix index for CADD indels |
| `revel` | `gs://YOUR_WORKSPACE_BUCKET/resources/revel/revel_grch38.tsv.gz` | REVEL GRCh38 tabix (prepare per Ensembl REVEL plugin docs) |
| `revel_index` | `gs://YOUR_WORKSPACE_BUCKET/resources/revel/revel_grch38.tsv.gz.tbi` | tabix index for REVEL |
| `dbnsfp` | `gs://YOUR_WORKSPACE_BUCKET/resources/dbnsfp/dbNSFP4.9c_grch38.gz` | VEP-ready dbNSFP (academic/commercial branch as licensed); supplies Eigen, MutationTaster, FATHMM |
| `dbnsfp_index` | `gs://YOUR_WORKSPACE_BUCKET/resources/dbnsfp/dbNSFP4.9c_grch38.gz.tbi` | tabix index for dbNSFP |

SIFT and PolyPhen also come from the VEP cache via `--sift b --polyphen b` (`vep_sift_polyphen`, default true) with no extra files. Override `dbnsfp_fields` if your dbNSFP build uses different column names (default includes `Eigen-phred_coding`, `MutationTaster_*`, `FATHMM_*`, and PolyPhen2/SIFT columns).

`LiftGnomadInfo` writes unified `INFO/gnomAD_AF` and `INFO/AF_popmax` as the **max** across genomes (`gnomADg_*`) and exomes (`gnomADe_*`). `pdvar` maps ClinVar `CLNSIG` to `clinvar_class`, and reads `CADD_PHRED` / `REVEL` / SIFT / PolyPhen / Eigen / MutationTaster / FATHMM from CSQ. Legacy `gnomad_vcf` still aliases genomes. Extra plugins (for example SpliceAI) can still use `vep_extra_args`.

## Optional inputs that are not in the JSON

These are unset in both JSON files. A default launch does not need them.

- **Introme GTF.** Stage GENCODE v44 to `gs://YOUR_WORKSPACE_BUCKET/resources/gencode/gencode.v44.annotation.gtf.gz` (see above). The Introme pin is `d6da48bc84abda66f588470f79abcf1244acfb21`.
- **Introme assets tar.** Optional `introme_assets_tar`. Unset, the fetch task downloads `https://github.com/CCICB/introme/archive/d6da48bc84abda66f588470f79abcf1244acfb21.tar.gz`.
- **VEP SpliceAI (and other plugins)** via `vep_extra_args` when you have the files. CADD and REVEL use the dedicated inputs above.
- **ExpansionHunter / Gauchian catalogs.** `eh_catalog` is already the vendored 19-locus hg38 catalog in this repo (ExpansionHunter 5.0.0). Gauchian 1.0.2 has no catalog file input.
