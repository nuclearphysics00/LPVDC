#!/bin/bash
set -euo pipefail

PBS_SCRIPT="worker_track_reco_v2.pbs"
#PBS_SCRIPT="worker_track_reco_optimize.pbs"


EXES=(
  "./track_reco_avalanche_pap_plate"
  "./track_reco_avalanche_pap_wirecath"
  "./track_reco_avalanche_ap_plate"
  "./track_reco_avalanche_ap_wirecath"
  "./track_reco_avalanche_a_plate"
  "./track_reco_avalanche_a_wirecath"
)

TAGS=(
  "pap_plate"
  "pap_wire"
  "ap_plate"
  "ap_wire"
  "a_plate"
  "a_wire"
)

if [ ! -f "$PBS_SCRIPT" ]; then
  echo "[ERROR] PBS script not found: $PBS_SCRIPT"
  exit 2
fi

if [ "${#EXES[@]}" -ne "${#TAGS[@]}" ]; then
  echo "[ERROR] EXES and TAGS size mismatch"
  exit 2
fi

for i in "${!EXES[@]}"; do
  exe="${EXES[$i]}"
  tag="${TAGS[$i]}"

  if [ ! -x "$exe" ]; then
    echo "[WARN] skip (not executable): $exe"
    continue
  fi

  echo "Submitting: tag=$tag  exe=$exe"
  qsub -v EXE="$exe",JOB_TAG="$tag" "$PBS_SCRIPT"
done