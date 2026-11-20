#!/usr/bin/env bash
# Build a panel-sliced ClinVar P/LP sites VCF.
#
# Re-run whenever the panel gene BED changes, then upload the output + .tbi to
# the Terra workspace bucket. No PdSample code changes are required for gene
# add/remove — only this extract and the staged File path.
#
# Example:
#   ./workflow/scripts/build_panel_clinvar.sh \
#     --bed resources/panel/monopdaus-genes.grch38.bed \
#     --clinvar gs://gcp-public-data--gnomad/resources/grch38/clinvar/clinvar_20250504.vcf.gz \
#     --out clinvar.20250504.panel.plp.vcf.bgz
#
# Public source used during dbNSFP work:
#   gs://gcp-public-data--gnomad/resources/grch38/clinvar/clinvar_20250504.vcf.gz
#
# Requires: bcftools, bgzip, tabix, and (for gs://) gsutil.

set -euo pipefail

BED=""
CLINVAR=""
OUT=""
WORKDIR=""
FLANK=0
RETRIES=4

usage() {
  sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
  echo "Options:"
  echo "  --bed PATH          Panel gene BED (chr-prefixed GRCh38)"
  echo "  --clinvar PATH      ClinVar sites VCF (.vcf.gz) or gs:// URI"
  echo "  --out PATH          Output .vcf.bgz path"
  echo "  --workdir PATH      Persistent work directory (default: mktemp, deleted on exit)"
  echo "  --flank INT         bp pad around each BED interval (default 0)"
  echo "  --retries INT       Attempts (default ${RETRIES})"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed) BED="$2"; shift 2 ;;
    --clinvar) CLINVAR="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --flank) FLANK="$2"; shift 2 ;;
    --retries) RETRIES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$BED" || -z "$CLINVAR" || -z "$OUT" ]]; then
  echo "ERROR: --bed, --clinvar, and --out are required" >&2
  usage >&2
  exit 2
fi
if [[ ! -s "$BED" ]]; then
  echo "ERROR: BED not found or empty: $BED" >&2
  exit 1
fi
command -v bcftools >/dev/null || { echo "ERROR: bcftools required" >&2; exit 1; }
command -v bgzip >/dev/null || { echo "ERROR: bgzip required" >&2; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix required" >&2; exit 1; }

KEEP_WORKDIR=0
if [[ -n "$WORKDIR" ]]; then
  KEEP_WORKDIR=1
  mkdir -p "$WORKDIR"
else
  WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/panel_clinvar.XXXXXX")"
fi
cleanup() {
  if [[ "$KEEP_WORKDIR" -eq 0 ]]; then
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

to_https() {
  local path="$1"
  if [[ "$path" == gs://* ]]; then
    echo "https://storage.googleapis.com/${path#gs://}"
  else
    echo "$path"
  fi
}

awk -v flank="$FLANK" 'BEGIN{OFS="\t"}
  $1 ~ /^#/ { next }
  NF >= 3 {
    chrom=$1; start=$2-flank; end=$3+flank;
    if (start < 0) start=0;
    if (end < start) end=start;
    print chrom, start, end
  }' "$BED" | sort -k1,1 -k2,2n > "$WORKDIR/regions.bed"

src_uri="$CLINVAR"
src="$(to_https "$src_uri")"
# Keep Pathogenic / Likely_pathogenic, including Pathogenic/Likely_pathogenic.
# Drop Conflicting_* and pure Benign/Likely_benign.
filter='(INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic") && INFO/CLNSIG!~"Conflicting" && INFO/CLNSIG!~"Benign"'

ok=0
attempt=1
while [[ "$attempt" -le "$RETRIES" ]]; do
  echo "Extracting ClinVar P/LP from $src (attempt ${attempt}/${RETRIES})" >&2
  rm -f "$WORKDIR/panel.plp.vcf.bgz" "$WORKDIR/panel.plp.vcf.bgz.tbi"
  set +e
  bcftools view -R "$WORKDIR/regions.bed" -i "$filter" -Oz -o "$WORKDIR/panel.plp.vcf.bgz" "$src" 2>"$WORKDIR/extract.err"
  status=$?
  set -e
  if [[ "$status" -eq 0 && -s "$WORKDIR/panel.plp.vcf.bgz" ]]; then
    tabix -f -p vcf "$WORKDIR/panel.plp.vcf.bgz"
    ok=1
    break
  fi
  echo "WARNING: stream extract failed (exit $status); trying local download" >&2
  sed -n '1,5p' "$WORKDIR/extract.err" >&2 || true
  if [[ "$src_uri" == gs://* ]] && command -v gsutil >/dev/null; then
    echo "Downloading $src_uri locally..." >&2
    set +e
    gsutil -m cp "$src_uri" "$src_uri.tbi" "$WORKDIR/" 2>"$WORKDIR/dl.err"
    dl_status=$?
    set -e
    local_src="$WORKDIR/$(basename "$src_uri")"
    if [[ "$dl_status" -eq 0 && -s "$local_src" ]]; then
      set +e
      bcftools view -R "$WORKDIR/regions.bed" -i "$filter" -Oz -o "$WORKDIR/panel.plp.vcf.bgz" "$local_src"
      status=$?
      set -e
      rm -f "$local_src" "${local_src}.tbi"
      if [[ "$status" -eq 0 && -s "$WORKDIR/panel.plp.vcf.bgz" ]]; then
        tabix -f -p vcf "$WORKDIR/panel.plp.vcf.bgz"
        ok=1
        break
      fi
    else
      echo "WARNING: gsutil download failed" >&2
      sed -n '1,5p' "$WORKDIR/dl.err" >&2 || true
    fi
  fi
  sleep $((attempt * 5))
  attempt=$((attempt + 1))
done

if [[ "$ok" -ne 1 ]]; then
  echo "ERROR: failed to extract ClinVar P/LP after ${RETRIES} attempts" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT")"
cp "$WORKDIR/panel.plp.vcf.bgz" "$OUT"
cp "$WORKDIR/panel.plp.vcf.bgz.tbi" "${OUT}.tbi"
echo "Wrote $OUT and ${OUT}.tbi" >&2
echo "Upload both to gs://YOUR_WORKSPACE_BUCKET/resources/clinvar/ and set PdSample.clinvar_plp_vcf / clinvar_plp_vcf_index." >&2
