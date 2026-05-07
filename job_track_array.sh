#!/bin/bash
#PBS -q AL
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -m abe
#PBS -M fumiya@rcnp.osaka-u.ac.jp
#PBS -N track_reco  
#PBS -t 1-100

# 安全運転
set -euo pipefail

# --- 配列ジョブのタスクID（シード値）を取得 ---
# PBS Pro / Torque の両方に対応できるようフォールバックを設ける
JOB_SEED=${PBS_ARRAYID:-${PBS_ARRAY_INDEX:-1}}

echo "start batch job: Task ID (Seed) = $JOB_SEED"

# --- 作業ディレクトリへ移動 ---
cd "$PBS_O_WORKDIR"
umask 0002

# --- 実行ファイル (qsub -v EXEC=... で上書き可能) ---
EXE=${EXEC:-./track_reco_avalanche_plate}

# --- T-L 曲線 CSV (qsub -v T_HIST_CSV=... で渡す、なければ自動検索) ---
if [[ -n "${T_HIST_CSV:-}" && -f "$T_HIST_CSV" ]]; then
  : # そのまま使う
elif [[ -n "${RUN_TAG:-}" ]]; then
  T_HIST_CSV="root/${RUN_TAG}/analysis_L0_prim/t_hist_nt.csv"
else
  # フォールバック: 最新 run の CSV を検索
  T_HIST_CSV=$(find root -name "t_hist_nt.csv" 2>/dev/null | sort | tail -1)
fi

if [[ -z "${T_HIST_CSV:-}" || ! -f "$T_HIST_CSV" ]]; then
  echo "[ERROR] t_hist_nt.csv が見つからない: ${T_HIST_CSV:-未指定}" >&2
  exit 1
fi

# --- ログセットアップ（親ジョブIDとタスクIDでディレクトリを分離） ---
BASE_JOBID=${PBS_JOBID%%.*}
BASE_JOBID=${BASE_JOBID%%\[*}
LOGDIR="logs/job_${BASE_JOBID}_task${JOB_SEED}"
mkdir -p "$LOGDIR"
: > "$LOGDIR/job.log"

# =================================================
# 環境設定（ROOT 6.32, GCC 10.3, Geant4, Garfield++）
# =================================================

# スタック拡張（自動配列のSEGV対策）
ulimit -s unlimited || true
ulimit -c 0 || true

# ROOT 6.32
source "/np1a/phanes/opt/root/6.32.06/bin/thisroot.sh"
PATH="$(echo "$PATH" | tr ':' '\n' | grep -v '/root_v6\.28\.06/' | paste -sd: -)"
LD_LIBRARY_PATH="$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v '/root_v6\.28\.06/' | paste -sd: -)"
export PATH LD_LIBRARY_PATH

# GCC 10.3
export LD_LIBRARY_PATH="/home/fumiya/local/gcc-10.3.0/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="/home/fumiya/local/gcc-10.3.0/bin${PATH:+:$PATH}"

# Geant4
source "$HOME/Geant4/geant4/install/bin/geant4.sh"

# Qt は使わない
unset QT_PLUGIN_PATH || true
unset QML2_IMPORT_PATH || true

# Garfield++
export GARFIELD_HOME="$HOME/garfieldpp-master/install"
export PATH="$GARFIELD_HOME${PATH:+:$PATH}"
export LD_LIBRARY_PATH="$GARFIELD_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# --- 環境情報をログ ---
{
  echo "[INFO] Job started at $(date)"
  echo "[INFO] Working dir: $PBS_O_WORKDIR"
  echo "[INFO] PWD=$(pwd)"
  [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ] && echo "[INFO] Node list: $(cat "$PBS_NODEFILE")"
  echo "[INFO] root-config --version: $(root-config --version 2>/dev/null || echo NA)"
} >> "$LOGDIR/job.log" 2>&1

# --- 実行ファイル存在確認 ---
if [ ! -x "$EXE" ]; then
  echo "[ERROR] Executable not found: $EXE" | tee -a "$LOGDIR/job.log"
  exit 2
fi

# =====================
# 実行＆エラーハンドリング
# =====================
# 第1引数に JOB_SEED、第2引数に T_HIST_CSV を渡して実行する
echo "[INFO] Running: $EXE $JOB_SEED $T_HIST_CSV" | tee -a "$LOGDIR/job.log"

set +e
"$EXE" "$JOB_SEED" "$T_HIST_CSV" >> "$LOGDIR/run_sim.out" 2>> "$LOGDIR/run_sim.err"
RUN_STATUS=$?
set -e

if [ $RUN_STATUS -eq 0 ]; then
  echo "[INFO] Job completed successfully at $(date)" >> "$LOGDIR/job.log"
  exit 0
fi

# 失敗時の詳細ログ
{
  echo "[ERROR] Exit code: $RUN_STATUS at $(date)"
  echo "----- tail run_sim.err -----"
  tail -n 80 "$LOGDIR/run_sim.err" || true
  echo "----- tail run_sim.out -----"
  tail -n 80 "$LOGDIR/run_sim.out" || true
} >> "$LOGDIR/job.log" 2>&1

exit $RUN_STATUS