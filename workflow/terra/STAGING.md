# Files to stage into the workspace bucket

`PdSample.example.json` and `PdSample.terra_inputs.json` point at public `gs://` objects wherever an anonymous copy exists. The inputs below do not. Their values stay `gs://YOUR_WORKSPACE_BUCKET/...`. Replace `YOUR_WORKSPACE_BUCKET` with the Terra workspace bucket when you launch. Do not commit a real bucket name, sample id, or billing project.

Mitochondria references, the GRCh38 FASTA, VerifyBamID 2.0.3 files, the vendored ExpansionHunter catalog, and the GATK-SV v0.29-beta resources are already public. They are not staged here. See [`docs/setup_checklist.md`](../../docs/setup_checklist.md).

## Placeholders

| Input | JSON | Placeholder | Upstream | Approx. size |
| --- | --- | --- | --- | --- |
| `PdSample.cram` | example only | `gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.cram` | your own GRCh38 CRAM | the alignment |
| `PdSample.crai` | example only | `gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.cram.crai` | index of that CRAM | small |
| `PdSample.vep_cache` | both | `gs://YOUR_WORKSPACE_BUCKET/resources/vep/homo_sapiens_vep_116_GRCh38.tar.gz` | Ensembl VEP cache release 116 | 26 GB (27,644,657,162 bytes) |
| `PdSample.somalier_sites` | both | `gs://YOUR_WORKSPACE_BUCKET/resources/somalier/sites.hg38.vcf.gz` | somalier v0.3.5 release notes, chr-prefixed hg38 | 260 KB (265,818 bytes) |

`PdSample.terra_inputs.json` binds the alignment to `${this.cram}` and `${this.crai}` on the sample table. The example JSON paths exist so Dockstore and `womtool` have a `File`. `SYNTHETIC_SAMPLE` is not a real sample.

`PdSample.requester_pays_project` stays `YOUR_GCP_PROJECT`. It is the GCP project that pays for requester-pays reads of `gs://gcp-public-data--broad-references` (the genome FASTA and the chrM bundle). It is not a file.

## Sample CRAM

Copy a CRAM aligned to `Homo_sapiens_assembly38` (chr-prefixed, dictionary MD5s matching that FASTA) and its `.crai` into the workspace bucket, then set the sample-table columns. This repository does not publish an example alignment.

```bash
gcloud storage cp /path/to/SYNTHETIC_SAMPLE.cram \
  gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.cram
gcloud storage cp /path/to/SYNTHETIC_SAMPLE.cram.crai \
  gs://YOUR_WORKSPACE_BUCKET/synthetic/SYNTHETIC_SAMPLE.cram.crai
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

Point `PdSample.vep_cache` and `PdSample.somalier_sites` at those two objects. Put the CRAM and CRAI on the sample table (`cram`, `crai`).

## Optional inputs that are not in the JSON

These are unset in both JSON files. A default launch does not need them.

- **Introme GTF.** `PdSample.introme_gtf` is already `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz` (GENCODE v44). The Introme pin is `d6da48bc84abda66f588470f79abcf1244acfb21`.
- **Introme assets tar.** Optional `introme_assets_tar`. Unset, the fetch task downloads `https://github.com/CCICB/introme/archive/d6da48bc84abda66f588470f79abcf1244acfb21.tar.gz`.
- **VEP plugin data** (CADD, REVEL, SpliceAI) passed through `vep_extra_args`. No public GCS copy was verified, and the files are license-restricted. Stage your own copies if you use them. They are not example-JSON inputs.
- **gnomAD and ClinVar** (`gnomad_vcf`, `clinvar_vcf`, and their indexes). Optional `--custom` annotations. They are not set in the example JSON.
- **ExpansionHunter / Gauchian catalogs.** `eh_catalog` is already the vendored 19-locus hg38 catalog in this repo (ExpansionHunter 5.0.0). Gauchian 1.0.2 has no catalog file input.
