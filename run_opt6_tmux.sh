#!/usr/bin/env bash
set -euo pipefail

# ===== 設定 =====
SESSION="opt6"

# 入力 CSV
INPUT="root/run_20260302_161622/analysis_L0_prim_exact/t_hist_nt.csv"

# ログ保存先
LOGDIR="logs_opt6_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOGDIR}"

# 実行する6バイナリ
JOBS=(
  "pap_plate     ./optimize_efield_pap_plate"
  "pap_wirecath  ./optimize_efield_pap_wirecath"
  "ap_plate      ./optimize_efield_ap_plate"
  "ap_wirecath   ./optimize_efield_ap_wirecath"
  "a_plate       ./optimize_efield_a_plate"
  "a_wirecath    ./optimize_efield_a_wirecath"
)

# 既存セッション確認
if tmux has-session -t "${SESSION}" 2>/dev/null; then
  echo "[Error] tmux session '${SESSION}' already exists."
  echo "Attach: tmux attach -t ${SESSION}"
  exit 1
fi

# バイナリ確認
for job in "${JOBS[@]}"; do
  name=$(echo "${job}" | awk '{print $1}')
  bin=$(echo "${job}" | awk '{print $2}')

  if [[ ! -x "${bin}" ]]; then
    echo "[Error] binary not found or not executable: ${bin}"
    echo "Try: chmod +x ${bin}"
    exit 1
  fi
done

# 最初の window を作成
first_name=$(echo "${JOBS[0]}" | awk '{print $1}')
tmux new-session -d -s "${SESSION}" -n "${first_name}"

# 各 window にコマンド投入
for i in "${!JOBS[@]}"; do
  name=$(echo "${JOBS[$i]}" | awk '{print $1}')
  bin=$(echo "${JOBS[$i]}" | awk '{print $2}')
  logfile="${LOGDIR}/${name}.log"

  if [[ "${i}" -ne 0 ]]; then
    tmux new-window -t "${SESSION}" -n "${name}"
  fi

  cmd="
echo '[Start] ${name}';
echo '[Binary] ${bin}';
echo '[Input] ${INPUT}';
echo '[Log] ${logfile}';
${bin} ${INPUT} > ${logfile} 2>&1;
status=\$?;
echo '[End] ${name}, status=' \$status;
echo 'Log: ${logfile}';
exec bash
"

  tmux send-keys -t "${SESSION}:${name}" "${cmd}" C-m
done

echo "[OK] Started tmux session: ${SESSION}"
echo "Attach:"
echo "  tmux attach -t ${SESSION}"
echo
echo "Log directory:"
echo "  ${LOGDIR}"
echo
echo "Watch logs:"
echo "  tail -f ${LOGDIR}/pap_plate.log"