#!/usr/bin/env bash
# Validate every WDL in this repo.
#
# Default (network required): womtool validates the committed files, including
# PdSample.wdl with its real HTTP imports:
#   gatk-sv v0.29-beta GATKSVPipelineSingleSample
#   GATK 4.7.0.0 MitochondriaPipeline
# miniwdl checks that same tree. It still substitutes stubs/MitochondriaPipeline.wdl
# because upstream AlignAndCall.wdl multiplies an optional Int? and miniwdl
# rejects that. GATK-SV is not stubbed. That mitochondria substitution is a
# miniwdl parser limit, not the offline fallback.
#
# --offline: no network. Both remote imports are replaced with the files in
# workflow/wdl/stubs/. Do not use this as the check that Dockstore or Terra run.
set -euo pipefail
offline=0
if [[ "${1:-}" == "--offline" ]]; then
  offline=1
elif [[ -n "${1:-}" ]]; then
  echo "usage: $0 [--offline]" >&2
  exit 2
fi

root=$(cd "$(dirname "$0")/../.." && pwd)
status=0

shopt -s nullglob
if [[ "$offline" -eq 0 ]]; then
  for wdl in "$root"/workflow/wdl/*.wdl "$root"/workflow/wdl/tasks/*.wdl "$root"/workflow/wdl/stubs/*.wdl; do
    echo "== womtool validate $wdl"
    if ! womtool validate "$wdl"; then
      status=1
    fi
  done
  work=$(mktemp -d)
  cp -a "$root/workflow/wdl/." "$work/"
  sed -i 's|https://raw.githubusercontent.com/broadinstitute/gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl|stubs/MitochondriaPipeline.wdl|' "$work/PdSample.wdl"
  for wdl in "$work"/*.wdl "$work"/tasks/*.wdl "$work"/stubs/*.wdl; do
    echo "== miniwdl check $wdl"
    if ! miniwdl check "$wdl"; then
      status=1
    fi
  done
  rm -rf "$work"
else
  work=$(mktemp -d)
  cp -a "$root/workflow/wdl/." "$work/"
  sed -i 's|https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.29-beta/wdl/GATKSVPipelineSingleSample.wdl|stubs/GATKSVPipelineSingleSample.wdl|' "$work/PdSample.wdl"
  sed -i 's|https://raw.githubusercontent.com/broadinstitute/gatk/4.7.0.0/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl|stubs/MitochondriaPipeline.wdl|' "$work/PdSample.wdl"
  for wdl in "$work"/*.wdl "$work"/tasks/*.wdl "$work"/stubs/*.wdl; do
    echo "== womtool validate $wdl"
    if ! womtool validate "$wdl"; then
      status=1
    fi
  done
  for wdl in "$work"/*.wdl "$work"/tasks/*.wdl "$work"/stubs/*.wdl; do
    echo "== miniwdl check $wdl"
    if ! miniwdl check "$wdl"; then
      status=1
    fi
  done
  rm -rf "$work"
fi
exit "$status"
