#!/bin/bash
set -euo pipefail

# ============================================================
# copy_analysis_macros_to_array_track.sh
#
# 目的:
#   カレントディレクトリにある解析マクロを、
#   root/run_*/vCat_*V/array_*/track/ 以下へコピーする。
#
# 使い方:
#   ./copy_analysis_macros_to_array_track.sh
#
# 特定の run だけにコピー:
#   ./copy_analysis_macros_to_array_track.sh root/run_20260322_205404_ap_wire
#
# 特定の array だけにコピー:
#   ./copy_analysis_macros_to_array_track.sh root/run_20260322_205404_ap_wire/vCat_-1500V/array_8401606
# ============================================================

# コピーしたいマクロ
MACROS=(
  "draw_diff_x.C"
  "draw_wire_ap_checks.C"
)

# 探索開始ディレクトリ
# 引数なしなら root 以下を探索
SEARCH_BASE="${1:-root}"

if [ ! -e "$SEARCH_BASE" ]; then
  echo "[ERROR] SEARCH_BASE not found: $SEARCH_BASE"
  exit 1
fi

# マクロが存在するか確認
missing=0
for macro in "${MACROS[@]}"; do
  if [ ! -f "$macro" ]; then
    echo "[ERROR] macro not found in current directory: $macro"
    missing=1
  fi
done

if [ "$missing" -ne 0 ]; then
  echo "[HINT] Put this script in the directory containing:"
  printf '  %s\n' "${MACROS[@]}"
  exit 1
fi

# コピー先 track ディレクトリを探す。
# 1) 引数が array_* なら、その直下の track を使う。
# 2) 引数が track なら、そのディレクトリを使う。
# 3) それ以外なら、配下の array_*/track を探す。
declare -a TRACK_DIRS=()

base_name="$(basename "$SEARCH_BASE")"

if [ "$base_name" = "track" ]; then
  TRACK_DIRS+=("$SEARCH_BASE")
elif [[ "$base_name" == array_* ]]; then
  mkdir -p "$SEARCH_BASE/track"
  TRACK_DIRS+=("$SEARCH_BASE/track")
else
  while IFS= read -r -d '' d; do
    TRACK_DIRS+=("$d")
  done < <(find "$SEARCH_BASE" -type d -path "*/array_*/track" -print0 | sort -z)
fi

if [ "${#TRACK_DIRS[@]}" -eq 0 ]; then
  echo "[WARN] no array_*/track directories found under: $SEARCH_BASE"
  echo "[HINT] If track does not exist yet, pass an array directory directly, e.g."
  echo "       ./copy_analysis_macros_to_array_track.sh root/run_xxx/vCat_-1500V/array_12345"
  exit 0
fi

echo "[INFO] copy macros:"
printf '  %s\n' "${MACROS[@]}"
echo "[INFO] number of destination track dirs = ${#TRACK_DIRS[@]}"

for track_dir in "${TRACK_DIRS[@]}"; do
  echo "[COPY] -> $track_dir"
  for macro in "${MACROS[@]}"; do
    cp -f "$macro" "$track_dir/"
  done
done

echo "[DONE] copied macros to array track directories."
