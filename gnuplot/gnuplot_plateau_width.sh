#!/bin/bash

#===========================================処理内容==============================================
# - ./output/に出力された各MHdata_xxx.csvの3行目と6行目からJ_red - J_greenの値とplateau widthの値を
#   抽出したファイルとしてplateau_width_dif_Jr_Jg.csvがある
# - 本スクリプトはplateau_width_dif_Jr_Jg.csv単一ファイルのデータをplotするためのものである
#=================================================================================================

#入力ファイル名
input_file1="../output/plateau_width_dif_Jg04_1of3.csv"
input_file2="../output/plateau_width_dif_Jg04_5of9.csv"
input_file3="../output/plateau_width_dif_Jg04_7of9.csv"

#出力ファイル名
output_filename="kagome_27site_dif_Jg04"

#======================入力ファイルが存在するかを確認する==============================
if [ ! -f "$input_file1" ]; then
    echo "FILE '$input_file1' does not exist. Please check the filename and try again."
fi

if [ ! -f "$input_file2" ]; then
    echo "FILE '$input_file2' does not exist. Please check the filename and try again."
fi

if [ ! -f "$input_file3" ]; then
    echo "FILE '$input_file3' does not exist. Please check the filename and try again."
fi

#=====================以下に入力ファイルが存在する場合の処理を書く=====================
if [ -n "$input_file1" ] && [ -n "$input_file2" ] && [ -n "$input_file3" ]; then
    echo "$input_file1 & $input_file2  & $input_file3 found."


    #同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
    if [ -f "${output_filename}.csv" ]; then
        echo "File already exists. Please enter a new filename({filename}.eps/png):"
        read new_filename
        output_filename="${new_filename}.eps"
    fi

    #グラフ名
    graph_name="Kagome 27siteplateau width J_g=0.4" 

    #gnuplotスクリプトファイルを作成
        gnuplot -p << EOF
        set terminal postscript eps color enhanced "Arial" 20
        set output "${output_filename}.eps"
        set multiplot layout 3,1
        set title "${graph_name}"
        set xlabel "J_r"
        set ylabel "plateau width"
        set key outside
        set title "1/3 plateau"
        plot "${input_file1}" using 1:2 with points title "M = 1/3"
        set title "5/9 plateau"
        plot "${input_file2}" using 1:2 with points title "M = 5/9"
        set title "7/9 plateau"
        plot "${input_file3}" using 1:2 with points title "M = 7/9"
        unset multiplot        

        set terminal pngcairo size 960, 960
        set output "${output_filename}.png"
        set multiplot layout 3,1 title "Kagome lattice 27 site"
        set xlabel "J_r"
        set ylabel "plateau width"
        set key outside
        set title "1/3 plateau"
        plot "${input_file1}" using 1:2 title "M = 1/3"
        set title "5/9 plateau"
        plot "${input_file2}" using 1:2 title "M = 5/9"
        set title "7/9 plateau"
        plot "${input_file3}" using 1:2 title "M = 7/9"
        unset multiplot
        set term qt 

EOF
else
    echo "$input_file1/2 not found."
    echo "Please specify an existing input file."
fi
