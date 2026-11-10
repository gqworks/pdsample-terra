#!/usr/bin/env bash
# Local WSL driver: build panel gnomAD AF into resources/gnomad/.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BED="$ROOT/resources/panel/monopdaus-genes.grch38.bed"
OUT_DIR="$ROOT/resources/gnomad"
OUT="$OUT_DIR/gnomad.joint.v4.1.panel.sites.vcf.bgz"
PARTS_DIR="/tmp/panel_gnomad_parts"
GNOMAD_DIR="gs://gcp-public-data--gnomad/release/4.1/vcf/joint"
LOG="$OUT_DIR/build_panel_gnomad_af.log"

mkdir -p "$OUT_DIR" "$PARTS_DIR"
chmod +x "$ROOT/workflow/scripts/build_panel_gnomad_af.sh"

{
  echo "ROOT=$ROOT"
  echo "BED=$BED"
  echo "OUT=$OUT"
  echo "PARTS_DIR=$PARTS_DIR"
  date -Is
  "$ROOT/workflow/scripts/build_panel_gnomad_af.sh" \
    --bed "$BED" \
    --gnomad-dir "$GNOMAD_DIR" \
    --out "$OUT" \
    --workdir "$PARTS_DIR" \
    --retries 5
  echo "DONE $(date -Is)"
  ls -lh "$OUT" "$OUT.tbi"
} 2>&1 | tee -a "$LOG"
