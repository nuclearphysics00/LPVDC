#!/bin/bash

# === ディレクトリ移動の記述を削除 ===
# 実行時のカレントディレクトリでそのまま処理する

# vP_の後の数値を基準にファイルをソートする
# sed でファイル名から数値を抽出し、数値としてソート(sort -n)して元のファイル名に戻す
# ※ -95, -90 ... -5, 0 のような「昇順」になります。
SORTED_FILES=$(ls *.png | sed -E 's/.*vP_(-?[0-9]+).*/\1 &/' | sort -k1,1n | cut -d' ' -f2-)

# 1. すべてのPNGを1つのPDFに結合する
echo "1/3: PNG 画像を vP の値でソートして PDF に結合中..."
convert $SORTED_FILES all_maps.pdf

# 2. 横向き・4枚割り付けのPDFを作成する
echo "2/3: 4 枚割り付け（横向き）PDF を作成中..."
pdfjam --landscape --nup 2x2 all_maps.pdf --outfile all_maps_4up_land.pdf

# 3. 両面印刷（短辺とじ）で出力する
echo "3/3: プリンター (A3C2570) へ送信中..."
lpr -P A3C2570 -o sides=two-sided-short-edge all_maps_4up_land.pdf

echo "印刷ジョブの送信が完了した。"