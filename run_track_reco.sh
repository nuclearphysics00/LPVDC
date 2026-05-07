#!/bin/bash
# run_track_reco.sh  — track_reco ステップのみを既存の run ディレクトリに対して実行するスクリプト
#
# 使い方:
#   zsh run_track_reco.sh [オプション]
#
# オプション:
#   -h  このヘルプを表示
#
# ============================================================
# CONFIGURATION
# 参照する run ディレクトリと実行パラメータをここで設定する
# ============================================================

## 参照する run ディレクトリ (root/ 以下)
## 例: root/run_20260322_202450  または絶対パスも可
RUN_DIR="root/run_20260322_202450"

## 使用する T-L 曲線 CSV
##   analysis_L0_prim : 1次電子のみ（通常はこちら）
##   analysis_L0_all  : 全アバランシェ電子
TL_SUBDIR="analysis_L0_prim"

## 検出器タイプ: pap | ap | a
DETECTOR="pap"

## ジオメトリ: wire | plate
GEOM="wire"

## 並列ジョブ数 (ローカル逐次実行なので「ジョブ ID」を 1 から N_JOBS まで順に回す)
N_JOBS=10

## 最初のジョブ ID (seed 兼 ID、連続番号で N_JOBS 本実行される)
START_SEED=1

# ============================================================
# CONFIGURATION エンド
# ============================================================

set -euo pipefail

# ヘルプ
if [[ "${1:-}" == "-h" ]]; then
  sed -n '2,10p' "$0"
  echo ""
  echo "CONFIGURATION はスクリプト冒頭の設定ブロックを直接編集すること。"
  exit 0
fi

# -----------------------------------------------------------------------
# 派生変数
# -----------------------------------------------------------------------
WORK="$(cd "$(dirname "$0")" && pwd)"
cd "$WORK"

export GARFIELD_HOME="${GARFIELD_HOME:-/Users/furukawafumiya/garfieldpp}"
export GARFIELD_INSTALL="${GARFIELD_HOME}/install"
export DYLD_LIBRARY_PATH="${GARFIELD_INSTALL}/lib:${DYLD_LIBRARY_PATH:-}"

# RUN_DIR を絶対パスに正規化
if [[ "$RUN_DIR" != /* ]]; then
  RUN_DIR="${WORK}/${RUN_DIR}"
fi

[[ "$GEOM" == "wire" ]] && GEOM_TAG="wirecath" || GEOM_TAG="plate"

TRACKRECO_BIN="${WORK}/track_reco_avalanche_${DETECTOR}_${GEOM_TAG}"
T_HIST_CSV="${RUN_DIR}/${TL_SUBDIR}/t_hist_nt.csv"

# -----------------------------------------------------------------------
# 事前チェック
# -----------------------------------------------------------------------
echo "============================================================"
echo " run_track_reco.sh  開始: $(date)"
echo "  RUN_DIR      = $RUN_DIR"
echo "  TL_SUBDIR    = $TL_SUBDIR"
echo "  DETECTOR     = $DETECTOR"
echo "  GEOM         = $GEOM_TAG"
echo "  N_JOBS       = $N_JOBS"
echo "  START_SEED   = $START_SEED"
echo "  BINARY       = $TRACKRECO_BIN"
echo "  T_HIST_CSV   = $T_HIST_CSV"
echo "============================================================"

if [[ ! -d "$RUN_DIR" ]]; then
  echo "[ERROR] RUN_DIR が見つからない: $RUN_DIR" >&2
  exit 1
fi

if [[ ! -f "$T_HIST_CSV" ]]; then
  echo "[ERROR] T-L 曲線 CSV が見つからない: $T_HIST_CSV" >&2
  echo "       先に analyze_drift.C を実行して CSV を生成すること" >&2
  exit 1
fi

if [[ ! -x "$TRACKRECO_BIN" ]]; then
  echo "[ERROR] バイナリが見つからない: $TRACKRECO_BIN" >&2
  echo "       compile_track_reco.sh を先に実行すること" >&2
  exit 1
fi

# -----------------------------------------------------------------------
# ログディレクトリ作成 (run ディレクトリの直下に配置)
# -----------------------------------------------------------------------
LOGDIR="${RUN_DIR}/logs_trackreco_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR"
echo "  LOGDIR       = $LOGDIR"
echo ""

# -----------------------------------------------------------------------
# track_reco を START_SEED から N_JOBS 本逐次実行
# -----------------------------------------------------------------------
echo "[track_reco] $N_JOBS ジョブを逐次実行中..."
FAIL=0
for (( i=0; i<N_JOBS; i++ )); do
  SEED=$(( START_SEED + i ))
  LOG="${LOGDIR}/track_job${SEED}.log"

  printf "  job %3d / %3d  (seed=%d) ..." "$((i+1))" "$N_JOBS" "$SEED"

  if "$TRACKRECO_BIN" "$SEED" "$T_HIST_CSV" > "$LOG" 2>&1; then
    echo " OK"
  else
    echo " FAIL (see $LOG)"
    FAIL=$(( FAIL + 1 ))
  fi
done

# -----------------------------------------------------------------------
# 完了メッセージ
# -----------------------------------------------------------------------
echo ""
echo "============================================================"
echo " track_reco 完了: $(date)"
echo "  成功: $(( N_JOBS - FAIL )) / $N_JOBS ジョブ"
if [[ $FAIL -gt 0 ]]; then
  echo "  失敗: $FAIL ジョブ  (ログ: $LOGDIR/)"
fi
echo ""
# 出力ファイルを確認
VCAT_DIR=$(ls -d "${RUN_DIR}"/vCat_*V 2>/dev/null | head -1 || true)
if [[ -n "$VCAT_DIR" ]]; then
  NROOT=$(find "${VCAT_DIR}/track" -name "track_results_job*.root" 2>/dev/null | wc -l | tr -d ' ')
  echo "  track ROOT ファイル数: $NROOT"
  echo "  出力先: ${VCAT_DIR}/track/"
fi
echo "  ログ: $LOGDIR/"
echo "============================================================"

[[ $FAIL -eq 0 ]] && exit 0 || exit 1
