#!/bin/bash
set -euo pipefail

OUT_ROOT="${1:-merged_track_results.root}"

shopt -s nullglob
files=( track_results_job*.root )
shopt -u nullglob

if [ "${#files[@]}" -eq 0 ]; then
  echo "[ERROR] track_results_job*.root が見つかりません"
  exit 1
fi

echo "[INFO] merging ${#files[@]} files -> $OUT_ROOT"
hadd -f "$OUT_ROOT" "${files[@]}"
echo "[DONE] wrote $OUT_ROOT"