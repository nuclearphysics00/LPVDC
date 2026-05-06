#!/bin/bash
# job_analyze.sh — Step 4: analyze_drift.C を PBS 上で実行する
#
# qsub -v で渡す環境変数:
#   RUN_TAG   : run_YYYYMMDD_HHMMSS  (必須)
#   OUTBASE   : root ディレクトリ名 (デフォルト: root)
#
# 依存: hadd (Step 3) の後に -W depend=afterok:<JID3> で投入すること
#
#PBS -q AL
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -N garfield_analyze
#PBS -m abe
#PBS -M fumiya@rcnp.osaka-u.ac.jp

set -euo pipefail

# --- 環境セットアップ ---
source /np1a/phanes/opt/root/6.32.06/bin/thisroot.sh
export PATH=$HOME/local/gcc-10.3.0/bin:$PATH
export LD_LIBRARY_PATH=$HOME/local/gcc-10.3.0/lib64:${LD_LIBRARY_PATH:-}
export PATH=$HOME/garfieldpp-master/install:$PATH
export GARFIELD_HOME=$HOME/garfieldpp-master/install
export LD_LIBRARY_PATH=$GARFIELD_HOME/lib64:${LD_LIBRARY_PATH:-}

cd "$PBS_O_WORKDIR"
umask 0002

# --- パラメータ ---
TAG=${RUN_TAG:?RUN_TAG is required (set via qsub -v)}
OUTBASE=${OUTBASE:-root}
OUTDIR="${OUTBASE}/${TAG}"
LOGDIR="logs/${TAG}"
MERGED="${OUTDIR}/grid_times.merged.root"

mkdir -p "$LOGDIR"

{
  echo "[INFO] analyze started at $(date)"
  echo "[INFO] PWD=$(pwd)  TAG=${TAG}  OUTDIR=${OUTDIR}"
} | tee -a "$LOGDIR/analyze.log"

# --- merged ROOT の存在確認 ---
if [[ ! -f "$MERGED" ]]; then
  echo "[ERROR] merged ROOT が見つかりません: $MERGED" | tee -a "$LOGDIR/analyze.log"
  exit 1
fi

# --- analyze_drift.C を探す ---
# 優先順位: (1) OUTDIR 内の analyze_drift.C, (2) 最新 run から検索
ANALYZE_C=""
if [[ -f "${OUTDIR}/analyze_drift.C" ]]; then
  ANALYZE_C="${OUTDIR}/analyze_drift.C"
else
  ANALYZE_C=$(find "${OUTBASE}" -name "analyze_drift.C" 2>/dev/null | sort | tail -1)
fi

if [[ -z "$ANALYZE_C" ]]; then
  echo "[ERROR] analyze_drift.C が見つかりません" | tee -a "$LOGDIR/analyze.log"
  exit 1
fi

# OUTDIR に analyze_drift.C がなければコピー
[[ "$(realpath "$ANALYZE_C")" != "$(realpath "${OUTDIR}/analyze_drift.C" 2>/dev/null || true)" ]] \
  && cp "$ANALYZE_C" "${OUTDIR}/analyze_drift.C"

echo "[INFO] analyze_drift.C = $ANALYZE_C" | tee -a "$LOGDIR/analyze.log"

# --- ROOT マクロ実行 ---
pushd "$OUTDIR" > /dev/null
root -l -b -q "analyze_drift.C" 2>&1 | tee -a "$PBS_O_WORKDIR/$LOGDIR/analyze.log"
popd > /dev/null

# --- 出力確認 ---
T_HIST="${OUTDIR}/analysis_L0_prim/t_hist_nt.csv"
if [[ ! -f "$T_HIST" ]]; then
  echo "[ERROR] t_hist_nt.csv が生成されませんでした: $T_HIST" | tee -a "$LOGDIR/analyze.log"
  exit 1
fi

echo "[INFO] analyze completed at $(date)  CSV=$T_HIST" | tee -a "$LOGDIR/analyze.log"
