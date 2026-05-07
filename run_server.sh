#!/bin/bash
# run_server.sh  — PBS バッチサーバー用パイプライン一括投入スクリプト
#
# 実行順序 (PBS 依存関係で自動直列化):
#   Step 1  optimize_efield_*   電圧最適化 [デフォルト実行 / -O でスキップ]
#             → 完了後に RMS_vs_Vp.csv を確認し、ベストの vPW 値を
#               main_fieldview_*.cc に手動で反映 → compile.sh を再実行
#               してから Step 2 を投入すること
#   Step 2  job_array.sh        fieldview アレイ (T-L 曲線元データ生成)
#   Step 3  job_merge.sh        hadd でシャードをマージ
#   Step 4  job_analyze.sh      analyze_drift.C で T-L 曲線 CSV 生成
#   Step 5  job_track_array.sh  track_reco アレイ (トラック再構成)
#
# 使い方:
#   bash run_server.sh [オプション]
#
# オプション:
#   -g GEOM       ジオメトリ: wire (デフォルト) または plate
#   -d DETECTOR   検出器タイプ: pap (デフォルト), ap, a
#   -n NEVENTS    シミュレーションイベント数 (デフォルト: 1000000)
#   -N N_SHARDS   fieldview 並列シャード数 (デフォルト: 100)
#   -T N_TRACK    track_reco 並列ジョブ数 (デフォルト: 100)
#   -S SEED_BASE  乱数シードのベース値 (デフォルト: 987654321)
#   -H HALO_MULT  ハロー判定係数 (デフォルト: 0)
#   -O            Step 1 (optimize_efield) をスキップする
#                 ※ すでに vPW を main_fieldview_*.cc に反映済みの場合に使用
#   -D            ドライラン: qsub を実行せずコマンドだけ表示する
#   -h            このヘルプを表示

set -euo pipefail

# -----------------------------------------------------------------------
# デフォルト値
# -----------------------------------------------------------------------
GEOM="wire"
DETECTOR="pap"
N_EVENTS=1000000
N_SHARDS=100
N_TRACK=100
SEED_BASE=987654321
HALO_MULT=0
RUN_OPTIMIZE=1   # デフォルトで Step 1 (optimize_efield) を実行する (-O でスキップ)
DRY_RUN=0

# -----------------------------------------------------------------------
# 引数パース
# -----------------------------------------------------------------------
while getopts "g:d:n:N:T:S:H:ODh" opt; do
  case $opt in
    g) GEOM="$OPTARG" ;;
    d) DETECTOR="$OPTARG" ;;
    n) N_EVENTS="$OPTARG" ;;
    N) N_SHARDS="$OPTARG" ;;
    T) N_TRACK="$OPTARG" ;;
    S) SEED_BASE="$OPTARG" ;;
    H) HALO_MULT="$OPTARG" ;;
    O) RUN_OPTIMIZE=0 ;;  # Step 1 をスキップしたい場合に指定
    D) DRY_RUN=1 ;;
    h)
      sed -n '2,25p' "$0"
      exit 0 ;;
    *) echo "[ERROR] 不明なオプション: -$OPTARG" >&2; exit 1 ;;
  esac
done

# -----------------------------------------------------------------------
# 派生変数
# -----------------------------------------------------------------------
WORK="$(cd "$(dirname "$0")" && pwd)"
cd "$WORK"

[[ "$GEOM" == "wire" ]] && GEOM_TAG="wirecath" || GEOM_TAG="plate"

ARRAY_END=$(( N_SHARDS - 1 ))

TAG="run_$(date +%Y%m%d_%H%M%S)"
OUTDIR="root/${TAG}"
LOGDIR="logs/${TAG}"
mkdir -p "$OUTDIR" "$LOGDIR"

FIELDVIEW_BIN="./fieldview_${DETECTOR}_${GEOM_TAG}"
OPTIMIZE_BIN="./optimize_efield_${DETECTOR}_${GEOM_TAG}"
TRACKRECO_BIN="./track_reco_avalanche_${DETECTOR}_${GEOM_TAG}"
T_HIST_CSV="${OUTDIR}/analysis_L0_prim/t_hist_nt.csv"

# qsub ラッパー (ドライランなら echo だけ)
qsub_run() {
  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[DRY] qsub $*" >&2
    echo "DRY_00000.hpc"
  else
    qsub "$@"
  fi
}

# -----------------------------------------------------------------------
# バイナリ存在確認
# -----------------------------------------------------------------------
MISSING=0
for bin in "$FIELDVIEW_BIN" "$TRACKRECO_BIN"; do
  [[ ! -x "$bin" ]] && echo "[ERROR] バイナリが見つからない: $bin" >&2 && MISSING=1
done
if [[ $MISSING -eq 1 ]]; then
  echo "       compile.sh / compile_track_reco.sh を先に実行すること" >&2
  exit 1
fi

echo "============================================================"
echo " run_server.sh  開始: $(date)"
echo "  DETECTOR  = $DETECTOR"
echo "  GEOM      = $GEOM_TAG"
echo "  N_EVENTS  = $N_EVENTS"
echo "  N_SHARDS  = $N_SHARDS"
echo "  N_TRACK   = $N_TRACK"
echo "  SEED_BASE = $SEED_BASE"
echo "  TAG       = $TAG"
[[ $DRY_RUN -eq 1 ]] && echo "  *** DRY RUN MODE ***"
echo "============================================================"

PREV_JID=""

# -----------------------------------------------------------------------
# Step 1: optimize_efield (デフォルト実行 / -O でスキップ)
# -----------------------------------------------------------------------
if [[ $RUN_OPTIMIZE -eq 1 ]]; then
  if [[ ! -x "$OPTIMIZE_BIN" ]]; then
    echo "[WARN] $OPTIMIZE_BIN が見つからない。Step 1 をスキップする。" >&2
  else
    echo ""
    echo "[Step 1] optimize_efield を投入中..."
    JID1=$(qsub_run \
      -N "opt_${TAG}" \
      -o "${LOGDIR}/opt.out" \
      -e "${LOGDIR}/opt.err" \
      -v "EXEC=${OPTIMIZE_BIN},RUN_TAG=${TAG},OUTBASE=root" \
      job_optimize.sh)
    echo "  [OK] JobID = ${JID1}"
    echo "  ※ optimize 完了後、RMS_vs_Vp.csv を確認して vPW を main_fieldview_*.cc に反映し"
    echo "     compile.sh を再実行してから Step 2 が自動投入されます。"
    PREV_JID="${JID1}"
  fi
else
  # Step 1 スキップ時 (-O 指定): すでに vPW をソースに反映済みの前提で Step 2 以降を投入
  # RUN_OPTIMIZE=0 で起動した場合はここを通る
  echo ""
  echo "[Step 1] optimize_efield はスキップ (-O が指定された、または RUN_OPTIMIZE=0 で起動)"
  echo "         vPW が main_fieldview_*.cc に正しく反映済みであることを確認すること。"
fi

# -----------------------------------------------------------------------
# Step 2: fieldview 配列ジョブ
# -----------------------------------------------------------------------
echo ""
echo "[Step 2] fieldview 配列ジョブを投入中 (${N_SHARDS} シャード)..."

DEPEND2=""
[[ -n "$PREV_JID" ]] && DEPEND2="-W depend=afterok:${PREV_JID}"

JID2=$(qsub_run \
  ${DEPEND2:+$DEPEND2} \
  -J "0-${ARRAY_END}" \
  -N "fv_${TAG}" \
  -o "${LOGDIR}/shard.\$PBS_ARRAY_INDEX.out" \
  -e "${LOGDIR}/shard.\$PBS_ARRAY_INDEX.err" \
  -v "N_SHARDS=${N_SHARDS},RUN_TAG=${TAG},EXEC=${FIELDVIEW_BIN},MAKE_MAPS=1,N_EVENTS=${N_EVENTS},SEED_BASE=${SEED_BASE},HALO_MULT=${HALO_MULT}" \
  job_array.sh)
echo "  [OK] JobID = ${JID2}"

# -----------------------------------------------------------------------
# Step 3: hadd マージ (配列全完了後)
# -----------------------------------------------------------------------
echo ""
echo "[Step 3] hadd マージジョブを投入中 (depend=afterokarray:${JID2})..."

JID3=$(qsub_run \
  -W "depend=afterokarray:${JID2}" \
  -N "merge_${TAG}" \
  -o "${LOGDIR}/merge.out" \
  -e "${LOGDIR}/merge.err" \
  -v "RUN_TAG=${TAG},OUTBASE=root" \
  job_merge.sh)
echo "  [OK] JobID = ${JID3}"

# -----------------------------------------------------------------------
# Step 4: analyze_drift.C (マージ完了後)
# -----------------------------------------------------------------------
echo ""
echo "[Step 4] analyze_drift ジョブを投入中 (depend=afterok:${JID3})..."

JID4=$(qsub_run \
  -W "depend=afterok:${JID3}" \
  -N "analyze_${TAG}" \
  -o "${LOGDIR}/analyze.out" \
  -e "${LOGDIR}/analyze.err" \
  -v "RUN_TAG=${TAG},OUTBASE=root" \
  job_analyze.sh)
echo "  [OK] JobID = ${JID4}"

# -----------------------------------------------------------------------
# Step 5: track_reco 配列ジョブ (CSV 生成完了後)
# -----------------------------------------------------------------------
echo ""
echo "[Step 5] track_reco 配列ジョブを投入中 (depend=afterok:${JID4})..."

JID5=$(qsub_run \
  -W "depend=afterok:${JID4}" \
  -J "1-${N_TRACK}" \
  -N "track_${TAG}" \
  -o "${LOGDIR}/track.\$PBS_ARRAY_INDEX.out" \
  -e "${LOGDIR}/track.\$PBS_ARRAY_INDEX.err" \
  -v "EXEC=${TRACKRECO_BIN},RUN_TAG=${TAG},OUTBASE=root,T_HIST_CSV=${T_HIST_CSV}" \
  job_track_array.sh)
echo "  [OK] JobID = ${JID5}"

# -----------------------------------------------------------------------
# サマリー
# -----------------------------------------------------------------------
echo ""
echo "============================================================"
echo " 全ジョブ投入完了: $(date)"
echo "  TAG              : ${TAG}"
echo "  OUTDIR           : ${OUTDIR}"
echo "  LOGDIR           : ${LOGDIR}"
if [[ $RUN_OPTIMIZE -eq 1 && -n "${JID1:-}" ]]; then
  echo "  Step1 optimize   : ${JID1}"
fi
echo "  Step2 fieldview  : ${JID2}"
echo "  Step3 merge      : ${JID3}"
echo "  Step4 analyze    : ${JID4}"
echo "  Step5 track_reco : ${JID5}"
echo ""
echo "  監視   : qstat -u \$(whoami)"
echo "  ログ   : ls ${LOGDIR}/"
echo "  結果   : ls ${OUTDIR}/"
echo "============================================================"

echo "${TAG}" > "${LOGDIR}/run_tag.txt"

