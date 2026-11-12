#!/usr/bin/env bash
# Validate every WDL in this repo.
#
# Default (network required for GATK-SV nested imports): womtool validates the
# committed files. PdSample imports a vendored gatk-sv v1.1.1
# GATKSVPipelineSingleSample (absolute nested HTTP imports) plus a vendored
# GATK 4.7.0.0 MitochondriaPipeline (fixed mt-scale disks; no remote mito fetch).
# miniwdl checks that same tree.
#
# --offline: no network. The vendored GATK-SV import is replaced with
# workflow/wdl/stubs/GATKSVPipelineSingleSample.wdl. Do not use this as the
# Dockstore/Terra check.
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
  wdl_globs=(
    "$root"/workflow/wdl/*.wdl
    "$root"/workflow/wdl/tasks/*.wdl
    "$root"/workflow/wdl/stubs/*.wdl
    "$root"/workflow/wdl/vendor/gatk-sv/*.wdl
    "$root"/workflow/wdl/vendor/mitochondria_m2_wdl/*.wdl
  )
  for wdl in "${wdl_globs[@]}"; do
    echo "== womtool validate $wdl"
    if ! womtool validate "$wdl"; then
      status=1
    fi
  done
  for wdl in "${wdl_globs[@]}"; do
    echo "== miniwdl check $wdl"
    if ! miniwdl check "$wdl"; then
      status=1
    fi
  done
else
  work=$(mktemp -d)
  cp -a "$root/workflow/wdl/." "$work/"
  sed -i 's|vendor/gatk-sv/GATKSVPipelineSingleSample.wdl|stubs/GATKSVPipelineSingleSample.wdl|' "$work/PdSample.wdl"
  offline_globs=(
    "$work"/*.wdl
    "$work"/tasks/*.wdl
    "$work"/stubs/*.wdl
    "$work"/vendor/gatk-sv/*.wdl
    "$work"/vendor/mitochondria_m2_wdl/*.wdl
  )
  for wdl in "${offline_globs[@]}"; do
    echo "== womtool validate $wdl"
    if ! womtool validate "$wdl"; then
      status=1
    fi
  done
  for wdl in "${offline_globs[@]}"; do
    echo "== miniwdl check $wdl"
    if ! miniwdl check "$wdl"; then
      status=1
    fi
  done
  rm -rf "$work"
fi
exit "$status"
