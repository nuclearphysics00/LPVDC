#!/bin/bash
set -eu

TAG=run_$(date +%Y%m%d_%H%M%S)
OUTDIR="root/${TAG}"
LOGDIR="logs/${TAG}"
mkdir -p "$OUTDIR" "$LOGDIR"

# 全シャード合計のイベント数（= N(t) 作るときの統計）
N_EVENTS=1000000   # 例えば 10 万イベント（必要に応じて変えてOK）

# 乱数シードのベース（再現したくなったらこれをメモ）
SEED_BASE=123456789

# ハロー判定用の係数（ワイヤ半径×係数）
HALO_MULT=0

qsub -J 0-99 \
  -N garf_plate_${TAG} \
  -o ${LOGDIR}/o.$TAG.\$PBS_ARRAY_INDEX \
  -e ${LOGDIR}/e.$TAG.\$PBS_ARRAY_INDEX \
  -v N_SHARDS=100,\
RUN_TAG=${TAG},\
EXEC=./fieldview_plate_threshold,\
MAKE_MAPS=1,\
N_EVENTS=${N_EVENTS},\
SEED_BASE=${SEED_BASE},\
HALO_MULT=${HALO_MULT} \
  job_array.sh
