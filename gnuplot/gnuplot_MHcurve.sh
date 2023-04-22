#!/bin/bash

#===========================================処理内容==============================================
# ./output/に出力されたM-h curveをplotするための複数のcsvファイルのデータを一つのepsファイルに出力する
# 判例はファイル名をそのまま利用するので、どういったパラメータに対するグラフ化は個人で覚えておく必要あり
# 変数output_file_nameは毎回変えないと上書きされてしまう点に注意
#=================================================================================================

#ファイル名のリストを取得する
file1="../output/Jg0/MHdata_param_red_11.csv"
file2="../output/Jg0/MHdata_param_red_26.csv"
file3="../output/Jg0/MHdata_param_red_35.csv"
file4="../output/Jg0/MHdata_param_red_50.csv"

#出力ファイル名
output_filename="./img_MH_curve/kagome_27site_Jg0.png"
#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
# if [ -f $ourput_filename ]; then
#     echo "File already exists. Please enter a new filename({filename.png}):"
#     read new_filename
#     output_filename="$new_filename"
# fi

#グラフ名
graph_name="kagome 27site J_g=0"

#gnuplotスクリプトファイルを作成
gnuplot -p << EOF
    set title "$graph_name"
    set xlabel "h"
    set ylabel "M/M_{sat}"
    set key outside top
    plot "$file1" using 1:2 with steps title "J_r = 0.1"
    replot "$file2" using 1:2 with steps dt (5,5) title "J_r = 0.5"
    replot "$file3" using 1:2 with steps dt (10,10) title "J_r = 0.68"
    replot "$file4" using 1:2 with steps dt (20,10) title "J_r = 0.98"
    set terminal pngcairo size 960, 960
    set output "$output_filename" ;replot
    set term qt
EOF