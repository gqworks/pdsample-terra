#!/usr/bin/env bash
# Fetch public resources whose URLs were checked. Does not download patient data.
# Somalier sites and the GATK-SV reference panel are not fetched: no verified URL.
set -euo pipefail

dest="${1:-workflow/resources}"
fetch_vep=0
if [[ "${2:-}" == "--vep" ]]; then
  fetch_vep=1
fi

mkdir -p "${dest}/verifybamid"

base="https://github.com/Griffan/VerifyBamID/raw/v2.0.3/resource"
prefix="1000g.phase3.100k.b38.vcf.gz.dat"
for suffix in UD V bed mu; do
  curl -fsSL -o "${dest}/verifybamid/${prefix}.${suffix}" "${base}/${prefix}.${suffix}"
done
echo "VerifyBamID files are ${dest}/verifybamid/${prefix}.{UD,V,bed,mu}."

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

Somalier sites: the v0.3.5 GitHub release has no sites VCF asset. Supply one
as somalier_sites. This script does not guess a URL.

GATK-SV v0.29-beta is required. Public hg38 and 1KG paths are already in
workflow/terra/PdSample.example.json (resolved from that tag's ref_panel_1kg.json
and resources_hg38.json). Terra attribute expressions are in
workflow/terra/PdSample.terra_inputs.json. This script does not download the panel.

Mitochondria references: the files MitochondriaPipeline 4.7.0.0 requires
(mt FASTA, shifted mt, BWA indexes, blacklist, chain, interval lists).
Stage the bundle you already use with that WDL.

Optional Introme assets tar and image mirrors: see docs/setup_checklist.md.
EOF
