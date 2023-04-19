#!/bin/bash

#===========================================処理内容==============================================
# ./output/に出力されたM-h curveをplotするための複数のcsvファイルのデータを一つのepsファイルに出力する
# 判例はファイル名をそのまま利用するので、どういったパラメータに対するグラフ化は個人で覚えておく必要あり
# 変数output_file_nameは毎回変えないと上書きされてしまう点に注意
#=================================================================================================

#ファイル名のリストを取得する
files=$(ls ./output/*.csv)

#出力ファイル名
output_filename="kaome_27site.eps"
#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
if [ -f $ourput_filename ]; then
    echo "File already exists. Please enter a new filename({filename}.eps):"
    read new_filename
    output_filename="${new_filename}.eps"
fi

#グラフ名
graph_name="kagome 27site"

#gnuplotスクリプトファイルを作成
gnuplot -p << EOF
    set term eps
    set output "$output_filename"
    set xlabel "h"
    set ylabel "M/M_{sat}"
    set key outside top
    plot for [file in "${files}"] file using 1:2 with linespoints title file
EOF