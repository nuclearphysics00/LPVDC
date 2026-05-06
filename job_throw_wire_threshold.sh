#!/bin/bash
set -eu

TAG=run_$(date +%Y%m%d_%H%M%S)
OUTDIR="root/${TAG}"
LOGDIR="logs/${TAG}"
mkdir -p "$OUTDIR" "$LOGDIR"

N_EVENTS=100000
SEED_BASE=987654321
HALO_MULT=0

qsub -J 0-99 \
  -N garf_wire_${TAG} \
  -o ${LOGDIR}/o.$TAG.\$PBS_ARRAY_INDEX \
  -e ${LOGDIR}/e.$TAG.\$PBS_ARRAY_INDEX \
  -v N_SHARDS=100,\
RUN_TAG=${TAG},\
EXEC=./fieldview_wirecath_threshold,\
MAKE_MAPS=1,\
N_EVENTS=${N_EVENTS},\
SEED_BASE=${SEED_BASE},\
HALO_MULT=${HALO_MULT} \
  job_array.sh
