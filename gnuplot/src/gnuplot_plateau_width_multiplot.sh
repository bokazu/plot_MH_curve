#!/bin/bash

#===========================================処理内容==============================================
# - ./output/に出力された各MHdata_xxx.csvの3行目と6行目からJ_red - J_greenの値とplateau widthの値を
#   抽出したファイルとしてplateau_width_dif_Jr_Jg.csvがある
# - 本スクリプトはJg00~Jg10のplateau_width_dif_Jr_Jg.csvのデータをmulti plotするためのものである
# - 1つのグラフに1/3、5/9、7/9プラトーのplateau幅の変化をまとめ、それを2行3列でまとめる
#=================================================================================================

#入力ファイルリスト名
J_numbers=(00 01 02 03 04 05 06 07 08 09 10)
for J_number in "${J_numbers[@]}"; do
    # echo $(ls ../../output/plateau_width/Jg${J_number}/)
    eval files${J_number}="($(ls -d ../../output/plateau_width/Jg${J_number}/*))"
done

echo "${files00[@]}"
echo "${files01[@]}"

#出力ファイル名
output_filename1="../img_plateau_width/kagome_27site_plateau_width_part1"
output_filename2="../img_plateau_width/kagome_27site_plateau_width_part2"

#=====================以下に入力ファイルが存在する場合の処理を書く=====================


    #同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
if [ -f "${output_filename1}.csv" ]; then
    echo "File already exists. Please enter a new filename(../img_plateau_width/{filename}.eps/png):"
    read new_filename
    output_filename1="../img_plateau_width/${new_filename}"
fi

if [ -f "${output_filename2}.csv" ]; then
    echo "File already exists. Please enter a new filename(../img_plateau_width/{filename}.eps/png):"
    read new_filename
    output_filename2="../img_plateau_width/${new_filename}"
fi

#グラフ名
graph_name1="Kagome 27siteplateau width part1"
graph_name2="Kagome 27siteplateau width part2" 

gnuplotスクリプトファイルを作成
    gnuplot -p << EOF
    load "string.plt"


    set terminal pngcairo size 1280, 960
    set output "${output_filename1}.png"
    set multiplot layout 2,3


    set title "${graph_name1}"
    set xlabel "J_r"
    set ylabel "plateau width"
    set key right top
    set size ratio -1

    itr=0
    set size 0.3,0.52
    set title "J_g=0.0"
    plot for [file in "${files00[@]}"] file using 1:2 with points title "1"

    set size 0.3,0.52
    set title "J_g=0.1"
    plot for [file in "${files01[@]}"] file using 1:2 with points title "2"

    set size 0.3,0.52
    set title "J_g=0.2"
    plot for [file in "${files02[@]}"] file using 1:2 with points title "3"
    
    set size 0.3,0.52
    set title "J_g=0.3"
    plot for [file in "${files03[@]}"] file using 1:2 with points title "4"

    set size 0.3,0.52
    set title "J_g=0.4"
    plot for [file in "${files04[@]}"] file using 1:2 with points title "5"

    set size 0.3,0.52
    set title "J_g=0.5"
    plot for [file in "${files05[@]}"] file using 1:2 with points title "5"

 
    unset multiplot
    set term qt 

    set terminal pngcairo size 1280, 960
    set output "${output_filename2}.png"
    set multiplot layout 2,3
    set title "${graph_name2}"
    set xlabel "J_r"
    set ylabel "plateau width"
    set key right top
    set size ratio -1

    set size 0.3,0.52
    set title "J_g=0.6"
    plot for [file in "${files06[@]}"] file using 1:2 with points title "1"

    set size 0.3,0.52
    set title "J_g=0.7"
    plot for [file in "${files07[@]}"] file using 1:2 with points title "2"

    set size 0.3,0.52
    set title "J_g=0.8"
    plot for [file in "${files08[@]}"] file using 1:2 with points title "3"
    
    set size 0.3,0.52
    set title "J_g=0.9"
    plot for [file in "${files09[@]}"] file using 1:2 with points title "4"

    set size 0.3,0.52
    set title "J_g=1.0"
    plot for [file in "${files10[@]}"] file using 1:2 with points title "5"
    unset multiplot
    set term qt

   
EOF
