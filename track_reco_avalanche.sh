#!/bin/bash
#PBS -q AL
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -m abe
#PBS -M fumiya@rcnp.osaka-u.ac.jp
#PBS -N track_reco  

# 安全運転
set -euo pipefail

echo "start batch job"

# --- 作業ディレクトリへ移動 ---
cd "$PBS_O_WORKDIR"
umask 0002

# --- ログセットアップ（JobID で衝突回避） ---
NOW=$(date +%Y%m%d_%H%M%S)
LOGDIR="logs/${NOW}_job${PBS_JOBID:-unknown}"
mkdir -p "$LOGDIR"
: > "$LOGDIR/job.log"

# --- 実行ファイル（今回実行するプログラムに変更） ---
EXE="./track_reco_avalanche_plate" # 実行ファイル

# =================================================
# 環境設定（ROOT 6.32 に統一、Qt 不使用、libstdc++ を明示）
# =================================================

# スタック拡張（自動配列のSEGV対策）、コア不要なら 0（必要なら unlimited）
ulimit -s unlimited || true
ulimit -c 0 || true

# ROOT 6.32 を明示して読み込む
source "/np1a/phanes/opt/root/6.32.06/bin/thisroot.sh"

# 旧 ROOT 6.28 パスを PATH/LD_LIBRARY_PATH から除去（念のため）
PATH="$(echo "$PATH" | tr ':' '\n' | grep -v '/root_v6\.28\.06/' | paste -sd: -)"
LD_LIBRARY_PATH="$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v '/root_v6\.28\.06/' | paste -sd: -)"
export PATH LD_LIBRARY_PATH

# GCC 10.3 を最優先（libstdc++ の取り違え防止）
export LD_LIBRARY_PATH="/home/fumiya/local/gcc-10.3.0/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="/home/fumiya/local/gcc-10.3.0/bin${PATH:+:$PATH}"

# Geant4（ビルド時と揃える）
source "$HOME/Geant4/geant4/install/bin/geant4.sh"

# Qt は使わない（混入防止）
unset QT_PLUGIN_PATH || true
unset QML2_IMPORT_PATH || true

# Garfield++
export GARFIELD_HOME="$HOME/garfieldpp-master/install"
export PATH="$GARFIELD_HOME${PATH:+:$PATH}"
export LD_LIBRARY_PATH="$GARFIELD_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# 参考（ビルド時に使ったなら情報として記録）
export CXXFLAGS="-march=native"
export CFLAGS="-march=native"

# --- 環境情報をログ ---
{
  echo "[INFO] Job started at $(date)"
  echo "[INFO] Working dir: $PBS_O_WORKDIR"
  echo "[INFO] PWD=$(pwd)"
  [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ] && echo "[INFO] Node list: $(cat "$PBS_NODEFILE")"
  echo "[INFO] ROOTSYS=$ROOTSYS"
  echo "[INFO] root-config --version: $(root-config --version 2>/dev/null || echo NA)"
  echo "[INFO] PATH=$PATH"
  echo "[INFO] LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
  echo "[INFO] Directory listing (top):"
  ls -l
} >> "$LOGDIR/job.log" 2>&1

# --- 実行ファイル存在確認 ---
if [ ! -x "$EXE" ]; then
  echo "[ERROR] Executable not found: $EXE" | tee -a "$LOGDIR/job.log"
  exit 2
fi

# --- 依存解決のチェックをログ保存 ---
{
  echo "[INFO] which g++: $(command -v g++ || echo 'not found')"
  echo "[INFO] ldd $EXE (Garfield/Geant4/ROOT parts)"
  ldd "$EXE" | egrep 'libGarfield|Geant4|libCore|libCling|stdc\+\+' || true
} >> "$LOGDIR/job.log" 2>&1

# =====================
# 実行＆エラーハンドリング
# =====================
echo "[INFO] Running: $EXE" | tee -a "$LOGDIR/job.log"

set +e
"$EXE" >> "$LOGDIR/run_sim.out" 2>> "$LOGDIR/run_sim.err"
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

# gdb があればバックトレース採取（-g ビルド推奨）
if command -v gdb >/dev/null 2>&1; then
  echo "[INFO] Collecting backtrace with gdb..." >> "$LOGDIR/job.log"
  gdb -q -batch -ex "set pagination off" -ex "run" -ex "bt" -ex "bt full" \
      -ex "info sharedlibrary" --args "$EXE" \
      >> "$LOGDIR/gdb_bt.txt" 2>&1 || true
fi

exit $RUN_STATUS