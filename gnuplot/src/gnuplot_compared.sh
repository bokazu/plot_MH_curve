#!/bin/bash

#===========================================処理内容===========================================================
# - ./output/に出力された各MHdata_xxx.csvの3行目と6行目からJ_red - J_redの値とplateau widthの値を
#   抽出したファイルとしてplateau_width_dif_Jr_Jg.csvがある
# - 本スクリプトはJg00~Jg10のplateau_width_dif_JgXX.csvとMHdata_param_XXX.csvのデータをmultiplotするためのものである
# - 1つのファイルの左側にJ_g=XXでのM=1/3,5/9,7/9プラトー幅の遷移の様子をplotしたグラフをplotする
# - 右側には上記のgraphにおいて何点かJ_r=??のデータを選びM-H curveをplotする 
#  - J_gの値ごとにM-H curveをplotしたいデータ点は異なるので、JgXXのディレクトリに対するloop処理は行わない 
#==============================================================================================================

J_green_num=07
J_green_val=0.7

#2. M-H curve plotのためのデータファイル
#plotしたいデータを選択 <-- MHdata_param_red_{J_red_num}.csv
length=5
J_red_num[1]=1  #J_red=0.0
J_red_num[2]=23 #J_red=0.6
J_red_num[3]=36 #J_red
J_red_num[4]=41 #J_red=0.8
J_red_num[5]=51 #J_red=1.0

#出力ファイル名
output_filename="../img_compared/Jg${J_green_num}/kagome_27site_compared_Jg${J_green_num}"

#入力ファイル名
#1. plateau width - J_redのグラフplotのためのデータファイル
eval width_files="($(ls ../../output/plateau_width/Jg${J_green_num}/*))"

echo "Selected data for plot plateau width are..."
for el in "${width_files[@]}"; do
    echo $el
done





echo "Selected data for plot MH curve are..."
for ((i=1; i <= length; i++)); do
    MH_files[i]="../../output/Jg${J_green_num}/MHdata_param_red_${J_red_num[${i}]}.csv"
    echo ${MH_files[$i]}
    J_red_val[i]=$(echo "scale=2; ${J_red_num[${i}]}*0.02-0.02" | bc | xargs printf "%.2f\n")
    echo ${J_red_val[$i]}
done

echo "output file name is..."
echo $output_filename

#グラフタイトル
graph_title1="Plateau width of 1/3, 5/9, 7/9 for various J_{red}"
graph_title2="M-H curve for each J_{red} values with J_g = $J_green_val"



#gnuplotスクリプト
gnuplot -p << EOF
load "string.plt"

#pngファイル
set terminal pngcairo size 1280, 1280
set output "${output_filename}.png"
set multiplot

set title "${graph_title2}"
set xlabel "h"
set ylabel "M/M_{sat}"
set key outside top title "J_{red}"
unset arrow
set xrange[0:3]
set yrange[0:1]
plot for [i=1:$length] word("${MH_files[@]}", i) using 1:2 with steps dt i lt i+13 lw 2 title word("${J_red_val[@]}",i)

set size 0.35, 0.45
set origin 0.1,0.5
set title "${graph_title1}"
set xlabel "J_r"
set ylabel "plateau width"
set key right top title "M/M_{sat}"
set xrange[0:1.8]
set yrange[0:1.6]
set size ratio -1


set arrow from ${J_red_val[1]},0 to ${J_red_val[1]},1.6 nohead dt(10,20) lt 14 lw 2
set arrow from ${J_red_val[2]},0 to ${J_red_val[2]},1.6 nohead dt(10,20) lt 15 lw 2
set arrow from ${J_red_val[3]},0 to ${J_red_val[3]},1.6 nohead dt(10,20) lt 16 lw 2
set arrow from ${J_red_val[4]},0 to ${J_red_val[4]},1.6 nohead dt(10,20) lt 17 lw 2
set arrow from ${J_red_val[5]},0 to ${J_red_val[5]},1.6 nohead dt(10,20) lt 18 lw 2
itr=-1
plot for [i=1:4] word("${width_files[@]}", i) using 1:2 with points pt i*3 title sprintf("%d/9", itr=itr+2)
# plot for [file in "${width_files[@]}"] file using 1:2 with point pt  title sprintf("%d/9", itr=itr+2)

unset multiplot
set term qt

EOF
