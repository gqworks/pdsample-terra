#!/usr/bin/env bash
# Fetch public resources whose URLs were checked. Does not download patient data.
# The GATK-SV reference panel is not fetched: the example JSON already points at
# the public objects. Mitochondria references are public GCS paths in that JSON.
set -euo pipefail

dest="${1:-workflow/resources}"
fetch_vep=0
if [[ "${2:-}" == "--vep" ]]; then
  fetch_vep=1
fi

mkdir -p "${dest}/verifybamid" "${dest}/somalier"

base="https://github.com/Griffan/VerifyBamID/raw/v2.0.3/resource"
prefix="1000g.phase3.100k.b38.vcf.gz.dat"
for suffix in UD V bed mu; do
  curl -fsSL -o "${dest}/verifybamid/${prefix}.${suffix}" "${base}/${prefix}.${suffix}"
done
echo "VerifyBamID files are ${dest}/verifybamid/${prefix}.{UD,V,bed,mu}."

# chr-prefixed hg38 sites linked from the somalier v0.3.5 release notes.
# sites.hg38.nochr.vcf.gz is the wrong build for this workflow.
curl -fL -o "${dest}/somalier/sites.hg38.vcf.gz" \
  "https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz"
echo "somalier sites VCF is ${dest}/somalier/sites.hg38.vcf.gz."
echo "Copy it to the workspace bucket. Commands are in workflow/terra/STAGING.md."

if [[ "${fetch_vep}" -eq 1 ]]; then
  mkdir -p "${dest}/vep"
  curl -fL -o "${dest}/vep/homo_sapiens_vep_116_GRCh38.tar.gz" \
    "https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz"
fi

cat <<'EOF'

Not downloaded here:

VEP cache (large; pass --vep to fetch it):
  https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz

Broad assembly38 (requester-pays; gsutil -u PROJECT):
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict

GENCODE v44 GTF (Introme):
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

ExpansionHunter catalog is vendored at workflow/assets/expansionhunter_pd.hg38.json
(19 loci copied from ExpansionHunter v5.0.0 hg38 variant_catalog.json).

GATK-SV v0.29-beta is required. Public hg38 and 1KG paths are already in
workflow/terra/PdSample.example.json (resolved from that tag's ref_panel_1kg.json
and resources_hg38.json). Terra attribute expressions are in
workflow/terra/PdSample.terra_inputs.json. This script does not download the panel.

Mitochondria references for MitochondriaPipeline 4.7.0.0 are already the public
gs://gcp-public-data--broad-references/hg38/v0/chrM/ objects in the example JSON
(same paths as that tag's ExampleInputsMitochondriaPipeline.json). This script
does not copy them.

Optional Introme assets tar and image mirrors: see docs/setup_checklist.md.
VEP cache and somalier staging commands: workflow/terra/STAGING.md.
EOF
