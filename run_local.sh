#!/bin/bash
# run_local.sh  — ローカル macOS 向けパイプライン一括実行スクリプト
#
# 実行順序:
#   Step 1  optimize_efield_*  (オプション: -o で有効化)
#   Step 2  fieldview_*        T-L 曲線の元データ生成
#   Step 3  hadd               シャードをマージ
#   Step 4  analyze_drift.C    T-L 曲線 CSV を生成
#   Step 5  track_reco_*       トラック再構成
#
# 使い方:
#   zsh run_local.sh [オプション]
#
# オプション:
#   -g GEOM       ジオメトリ: wire (デフォルト) または plate
#   -d DETECTOR   検出器タイプ: pap (デフォルト), ap, a
#   -n NEVENTS    シミュレーションイベント数 (デフォルト: 50000)
#   -j JOB_SEED   track_reco のジョブシード (デフォルト: 0)
#   -o            Step 1 の optimize_efield も実行する
#   -h            このヘルプを表示

set -euo pipefail

# -----------------------------------------------------------------------
# デフォルト値
# -----------------------------------------------------------------------
GEOM="wire"          # wire | plate
DETECTOR="pap"       # pap | ap | a
N_EVENTS=50000
JOB_SEED=0
RUN_OPTIMIZE=0

# -----------------------------------------------------------------------
# 引数パース
# -----------------------------------------------------------------------
while getopts "g:d:n:j:oh" opt; do
  case $opt in
    g) GEOM="$OPTARG" ;;
    d) DETECTOR="$OPTARG" ;;
    n) N_EVENTS="$OPTARG" ;;
    j) JOB_SEED="$OPTARG" ;;
    o) RUN_OPTIMIZE=1 ;;
    h)
      sed -n '2,20p' "$0"
      exit 0 ;;
    *) echo "[ERROR] 不明なオプション: -$OPTARG" >&2; exit 1 ;;
  esac
done

# -----------------------------------------------------------------------
# 派生変数
# -----------------------------------------------------------------------
WORK="$(cd "$(dirname "$0")" && pwd)"
cd "$WORK"

export GARFIELD_HOME="${GARFIELD_HOME:-/Users/furukawafumiya/garfieldpp}"
export GARFIELD_INSTALL="${GARFIELD_HOME}/install"
export DYLD_LIBRARY_PATH="${GARFIELD_INSTALL}/lib:${DYLD_LIBRARY_PATH:-}"

if [[ "$GEOM" == "wire" ]]; then
  GEOM_TAG="wirecath"
  GEOM_DEF="GEOM_WIRE"
else
  GEOM_TAG="plate"
  GEOM_DEF="GEOM_PLATE"
fi

TAG="run_$(date +%Y%m%d_%H%M%S)"
OUTDIR="root/${TAG}"
LOGDIR="logs/${TAG}"
mkdir -p "$OUTDIR" "$LOGDIR"

FIELDVIEW_BIN="./fieldview_${DETECTOR}_${GEOM_TAG}"
OPTIMIZE_BIN="./optimize_efield_${DETECTOR}_${GEOM_TAG}"
TRACKRECO_BIN="./track_reco_avalanche_${DETECTOR}_${GEOM_TAG}"

echo "============================================================"
echo " run_local.sh  開始: $(date)"
echo "  DETECTOR  = $DETECTOR"
echo "  GEOM      = $GEOM_TAG"
echo "  N_EVENTS  = $N_EVENTS"
echo "  JOB_SEED  = $JOB_SEED"
echo "  OUTDIR    = $OUTDIR"
echo "============================================================"

# -----------------------------------------------------------------------
# Step 0: ビルド確認
# -----------------------------------------------------------------------
echo ""
echo "[Step 0] バイナリの存在確認..."
MISSING=0
for bin in "$FIELDVIEW_BIN" "$TRACKRECO_BIN"; do
  if [[ ! -x "$bin" ]]; then
    echo "  [WARN] バイナリが見つかりません: $bin"
    echo "         compile.sh / compile_track_reco.sh を先に実行してください。"
    MISSING=1
  fi
done
if [[ $MISSING -eq 1 ]]; then
  exit 1
fi
echo "  [OK] バイナリ確認完了"

# -----------------------------------------------------------------------
# Step 1 (オプション): optimize_efield
# -----------------------------------------------------------------------
if [[ $RUN_OPTIMIZE -eq 1 ]]; then
  echo ""
  echo "[Step 1] optimize_efield 実行中..."
  if [[ ! -x "$OPTIMIZE_BIN" ]]; then
    echo "  [WARN] $OPTIMIZE_BIN が見つかりません。compile_optimize_efield.sh を実行してください。"
  else
    "$OPTIMIZE_BIN" 2>&1 | tee "$LOGDIR/optimize.out"
    echo "  [OK] optimize_efield 完了"
  fi
else
  echo ""
  echo "[Step 1] optimize_efield はスキップします (-o で有効化)"
fi

# -----------------------------------------------------------------------
# Step 2: fieldview でドリフトシミュレーション → grid_times.shard0000.root
# -----------------------------------------------------------------------
echo ""
echo "[Step 2] fieldview ドリフトシミュレーション実行中..."
export N_SHARDS=1
export SHARD_ID=0
export MAKE_MAPS=1
export N_EVENTS
export OUT_ROOT="${OUTDIR}/grid_times.shard0000.root"
export HALO_MULT=0

"$FIELDVIEW_BIN" 2>&1 | tee "$LOGDIR/fieldview.out"
echo "  [OK] fieldview 完了 → $OUT_ROOT"

# -----------------------------------------------------------------------
# Step 3: hadd でシャードをマージ
# -----------------------------------------------------------------------
echo ""
echo "[Step 3] ROOT ファイルをマージ中..."
MERGED="${OUTDIR}/grid_times.merged.root"
hadd -f "$MERGED" "${OUTDIR}"/grid_times.shard*.root
echo "  [OK] マージ完了 → $MERGED"

# -----------------------------------------------------------------------
# Step 4: analyze_drift.C で T-L 曲線 CSV を生成
# -----------------------------------------------------------------------
echo ""
echo "[Step 4] analyze_drift.C で T-L 曲線 CSV を生成中..."
ANALYZE_SRC="${WORK}/root/run_20260322_141320/analyze_drift.C"
if [[ ! -f "$ANALYZE_SRC" ]]; then
  # フォールバック: 最新の run から探す
  ANALYZE_SRC=$(find "${WORK}/root" -name "analyze_drift.C" | sort | tail -1)
fi
if [[ -z "$ANALYZE_SRC" ]]; then
  echo "  [ERROR] analyze_drift.C が見つかりません"
  exit 1
fi
cp "$ANALYZE_SRC" "$OUTDIR/analyze_drift.C"

pushd "$OUTDIR" > /dev/null
root -l -b -q "analyze_drift.C" 2>&1 | tee "$WORK/$LOGDIR/analyze_drift.out"
popd > /dev/null

T_HIST_CSV="${OUTDIR}/analysis_L0_prim/t_hist_nt.csv"
if [[ ! -f "$T_HIST_CSV" ]]; then
  echo "  [ERROR] t_hist_nt.csv が生成されませんでした"
  exit 1
fi
echo "  [OK] T-L 曲線 CSV 完了 → $T_HIST_CSV"

# -----------------------------------------------------------------------
# Step 5: track_reco でトラック再構成
# -----------------------------------------------------------------------
echo ""
echo "[Step 5] track_reco_avalanche 実行中..."
"$TRACKRECO_BIN" "$JOB_SEED" "$T_HIST_CSV" 2>&1 | tee "$LOGDIR/track_reco.out"
echo "  [OK] track_reco 完了"

# -----------------------------------------------------------------------
# 完了メッセージ
# -----------------------------------------------------------------------
echo ""
echo "============================================================"
echo " 全 Step 完了: $(date)"
echo "  ログ     : $LOGDIR/"
echo "  ROOT出力 : $OUTDIR/"
echo "  T-L CSV  : $T_HIST_CSV"
echo "============================================================"
