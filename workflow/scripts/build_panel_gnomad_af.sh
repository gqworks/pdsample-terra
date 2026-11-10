#!/usr/bin/env bash
# Build a panel-sliced gnomAD v4.1 joint sites VCF (AF_joint + AF_grpmax_joint).
#
# Re-run whenever the panel gene BED changes, then upload the output + .tbi to
# the Terra workspace bucket. No PdSample code changes are required for gene
# add/remove — only this extract and the staged File path.
#
# Example:
#   ./workflow/scripts/build_panel_gnomad_af.sh \
#     --bed resources/panel/monopdaus-genes.grch38.bed \
#     --gnomad-dir gs://gcp-public-data--gnomad/release/4.1/vcf/joint \
#     --out gnomad.joint.v4.1.panel.sites.vcf.bgz \
#     --workdir /tmp/panel_gnomad_parts
#
# Public source (per chromosome):
#   gs://gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chrN.vcf.bgz
#
# Requires: bcftools, bgzip, tabix, and (for gs://) gsutil.

set -euo pipefail

BED=""
GNOMAD_DIR=""
OUT=""
WORKDIR=""
FLANK=0
RETRIES=4
PREFIX="gnomad.joint.v4.1.sites"
KEEP_INFO="^INFO/AF_joint,INFO/AF_grpmax_joint,INFO/AC_joint,INFO/AN_joint,INFO/nhomalt_joint"

usage() {
  sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
  echo "Options:"
  echo "  --bed PATH          Panel gene BED (chr-prefixed GRCh38)"
  echo "  --gnomad-dir PATH   Directory or gs:// prefix with ${PREFIX}.chr*.vcf.bgz"
  echo "  --out PATH          Output .vcf.bgz path"
  echo "  --workdir PATH      Persistent parts directory (default: mktemp, deleted on exit)"
  echo "  --flank INT         bp pad around each BED interval (default 0)"
  echo "  --retries INT       Attempts per chromosome (default ${RETRIES})"
  echo "  --prefix STR        Chrom file prefix (default ${PREFIX})"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed) BED="$2"; shift 2 ;;
    --gnomad-dir) GNOMAD_DIR="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --flank) FLANK="$2"; shift 2 ;;
    --retries) RETRIES="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$BED" || -z "$GNOMAD_DIR" || -z "$OUT" ]]; then
  echo "ERROR: --bed, --gnomad-dir, and --out are required" >&2
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
  WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/panel_gnomad.XXXXXX")"
fi
cleanup() {
  if [[ "$KEEP_WORKDIR" -eq 0 ]]; then
    rm -rf "$WORKDIR"
  fi
}
trap cleanup EXIT

# Prefer HTTPS for public GCS objects (more retry-friendly than gs:// via htslib).
to_https() {
  local path="$1"
  if [[ "$path" == gs://* ]]; then
    echo "https://storage.googleapis.com/${path#gs://}"
  else
    echo "$path"
  fi
}

# Normalize BED to chrom / start / end; optional flank; sort for bcftools -R.
awk -v flank="$FLANK" 'BEGIN{OFS="\t"}
  $1 ~ /^#/ { next }
  NF >= 3 {
    chrom=$1; start=$2-flank; end=$3+flank;
    if (start < 0) start=0;
    if (end < start) end=start;
    print chrom, start, end
  }' "$BED" | sort -k1,1 -k2,2n > "$WORKDIR/regions.bed"

mapfile -t CHROMS < <(cut -f1 "$WORKDIR/regions.bed" | sed 's/^chr//' | sort -u)
PARTS=()
for bare in "${CHROMS[@]}"; do
  [[ -z "$bare" ]] && continue
  chrom="chr${bare}"
  src_uri="${GNOMAD_DIR%/}/${PREFIX}.${chrom}.vcf.bgz"
  src="$(to_https "$src_uri")"
  awk -v c="$chrom" -v b="chr${bare}" '($1==c || $1==b || $1=="'"$bare"'")' \
    "$WORKDIR/regions.bed" > "$WORKDIR/${chrom}.bed"
  if [[ ! -s "$WORKDIR/${chrom}.bed" ]]; then
    continue
  fi
  part="$WORKDIR/${chrom}.slim.vcf.bgz"
  if [[ -s "$part" && -s "${part}.tbi" ]]; then
    echo "Reusing existing $part" >&2
    PARTS+=("$part")
    continue
  fi

  ok=0
  attempt=1
  while [[ "$attempt" -le "$RETRIES" ]]; do
    echo "Extracting $chrom from $src (attempt ${attempt}/${RETRIES})" >&2
    rm -f "$part" "${part}.tbi" "$WORKDIR/${chrom}.full.vcf.bgz" "$WORKDIR/${chrom}.full.vcf.bgz.tbi"
    set +e
    bcftools view -R "$WORKDIR/${chrom}.bed" -Ou "$src" 2>"$WORKDIR/${chrom}.err" \
      | bcftools annotate -x "$KEEP_INFO" -Oz -o "$part"
    status=$?
    set -e
    if [[ "$status" -eq 0 && -s "$part" ]]; then
      tabix -f -p vcf "$part"
      ok=1
      break
    fi
    echo "WARNING: stream extract failed for $chrom (exit $status); trying local download" >&2
    if [[ -s "$WORKDIR/${chrom}.err" ]]; then
      sed -n '1,5p' "$WORKDIR/${chrom}.err" >&2 || true
    fi
    rm -f "$part" "${part}.tbi"

    # Fallback: copy the chrom VCF locally (needs ~70 GiB free), then extract.
    if [[ "$src_uri" == gs://* ]] && command -v gsutil >/dev/null; then
      echo "Downloading $src_uri locally..." >&2
      set +e
      gsutil -m cp "$src_uri" "$src_uri.tbi" "$WORKDIR/" 2>"$WORKDIR/${chrom}.dl.err"
      dl_status=$?
      set -e
      local_src="$WORKDIR/$(basename "$src_uri")"
      if [[ "$dl_status" -eq 0 && -s "$local_src" ]]; then
        set +e
        bcftools view -R "$WORKDIR/${chrom}.bed" -Ou "$local_src" \
          | bcftools annotate -x "$KEEP_INFO" -Oz -o "$part"
        status=$?
        set -e
        rm -f "$local_src" "${local_src}.tbi"
        if [[ "$status" -eq 0 && -s "$part" ]]; then
          tabix -f -p vcf "$part"
          ok=1
          break
        fi
      else
        echo "WARNING: gsutil download failed for $chrom" >&2
        sed -n '1,5p' "$WORKDIR/${chrom}.dl.err" >&2 || true
      fi
    fi

    sleep $((attempt * 5))
    attempt=$((attempt + 1))
  done

  if [[ "$ok" -ne 1 ]]; then
    echo "ERROR: failed to extract $chrom after ${RETRIES} attempts" >&2
    exit 1
  fi
  PARTS+=("$part")
done

if [[ ${#PARTS[@]} -eq 0 ]]; then
  echo "ERROR: no chromosome extracts produced" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT")"
if [[ ${#PARTS[@]} -eq 1 ]]; then
  cp "${PARTS[0]}" "$OUT"
  cp "${PARTS[0]}.tbi" "${OUT}.tbi"
else
  bcftools concat -a -Oz -o "$OUT" "${PARTS[@]}"
  tabix -f -p vcf "$OUT"
fi

echo "Wrote $OUT and ${OUT}.tbi" >&2
echo "Upload both to gs://YOUR_WORKSPACE_BUCKET/resources/gnomad/ and set PdSample.gnomad_af_vcf / gnomad_af_vcf_index." >&2
