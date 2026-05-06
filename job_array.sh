#!/bin/bash
#PBS -q AL
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -N garfield_shard
#PBS -m abe
#PBS -M fumiya@rcnp.osaka-u.ac.jp

set -euo pipefail

# --- 環境セットアップ ---
source /np1a/phanes/opt/root/6.32.06/bin/thisroot.sh
export PATH=$HOME/garfieldpp-master/install:$PATH
export GARFIELD_HOME=$HOME/garfieldpp-master/install
export LD_LIBRARY_PATH=$GARFIELD_HOME/lib64:$LD_LIBRARY_PATH

cd "$PBS_O_WORKDIR"
umask 0002

# --- 環境変数（qsub -v から来る想定） ---
EXEC=${EXEC:-./fieldview_wirecath}   # plate のときは外から EXEC=./fieldview_plate
N_SHARDS=${N_SHARDS:-100}
RUN_TAG=${RUN_TAG:-run_$(date +%Y%m%d_%H%M%S)}
SHARD_ID=${PBS_ARRAY_INDEX:?ARRAY only}

OUTDIR="root/${RUN_TAG}"
LOGDIR="logs/${RUN_TAG}"
mkdir -p "$OUTDIR" "$LOGDIR"

# 電場マップは shard 0 のときだけ
if [ "$SHARD_ID" -eq 0 ]; then
  MAKE_MAPS=${MAKE_MAPS:-1}
else
  MAKE_MAPS=0
fi

# ここで N_EVENTS / SEED_BASE / HALO_MULT はそのまま渡すだけ
# （指定がなければ実行バイナリ側のデフォルトを使う）
export N_SHARDS SHARD_ID MAKE_MAPS N_EVENTS SEED_BASE HALO_MULT

# 出力 ROOT ファイル名
export OUT_ROOT="$OUTDIR/grid_times.shard$(printf "%04d" "$SHARD_ID").root"

echo "[INFO] shard=$SHARD_ID N_SHARDS=$N_SHARDS OUT=$OUT_ROOT EXEC=$EXEC"
echo "[INFO] N_EVENTS=${N_EVENTS:-default} SEED_BASE=${SEED_BASE:-default} HALO_MULT=${HALO_MULT:-3.0}"

"$EXEC" >> "$LOGDIR/shard_$(printf "%04d" "$SHARD_ID").out" \
        2>> "$LOGDIR/shard_$(printf "%04d" "$SHARD_ID").err"
