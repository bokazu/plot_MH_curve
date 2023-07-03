#!/bin/bash

#===========================================処理内容==============================================
# ./output/に出力されたM-h curveをplotするための複数のcsvファイルのデータを一つのepsファイルに出力する
# 判例はファイル名をそのまま利用するので、どういったパラメータに対するグラフ化は個人で覚えておく必要あり
# 変数output_file_nameは毎回変えないと上書きされてしまう点に注意
#=================================================================================================

#ファイル名のリストを取得する
# J_num="Jg00"
# J_val="J_g=0.0"


file1="../../output/kagome/27site/output_param1/MHdata.csv"
# file2="../../output/kagome/27site/trimer/MHdata.csv"

# file1="../output/${J_num}/MHdata_param_red_11.csv"
# file2="../output/${J_num}/MHdata_param_red_26.csv"
# file3="../output/${J_num}/MHdata_param_red_35.csv"
# file4="../output/${J_num}/MHdata_param_red_50.csv"

#出力ファイル名
output_filename="../img_MH_curve/kagome_27site_hexagon_and_trimer"
# output_filename="./img_MH_curve/${J_num}/kagome_27site_${J_num}.png"
#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
# if [ -f $ourput_filename ]; then
#     echo "File already exists. Please enter a new filename({filename.png}):"
#     read new_filename
#     output_filename="$new_filename"
# fi

#グラフ名
graph_name="27 site, hexagon"
# graph_name="kagome 27site ${J_val}"

#gnuplotスクリプトファイルを作成
gnuplot -p << EOF
    set terminal pngcairo size 1280, 960
    set output "${output_filename}.png"

    set title "${graph_name}"
    set xlabel "h"
    set ylabel "M/M_{sat}"
    set key outside top
    set xrange [0:3.5]
    set yrange [0:1]

    plot "${file1}" using 1:2 with steps lw 2 linetype rgbcolor 'blue' title "hexagon","${file2}" using 1:2 with steps dt(15,20) lw 2 linetype rgbcolor 'red' title "trimer"
    # ,"${file3}" using 1:2 with steps dt(15,20) lw 2 linetype rgbcolor 'red' title "36site(J_r = 0.86, J_g = 0.06. J_b=1.0)"
    
    # replot "$file2" using 1:2 with steps dt (5,5) title "J_r = 0.5"
    # replot "$file3" using 1:2 with steps dt (10,10) title "J_r = 0.68"
    # replot "$file4" using 1:2 with steps dt (20,10) title "J_r = 0.98"
    set term qt
EOF
