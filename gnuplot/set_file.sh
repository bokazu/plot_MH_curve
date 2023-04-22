#!/bin/bash


#===================相互作用とplateau幅の関係を調べるための処理===
plateau_width_output_file1="../output/plateau_width_dif_Jg04_1of3.csv" #[コメント]projectタイトルのような変数を用意して、それを末尾にくっつけると良いかも?
plateau_width_output_file2="../output/plateau_width_dif_Jg04_5of9.csv"
plateau_width_output_file3="../output/plateau_width_dif_Jg04_7of9.csv"

#=============================1/3プラトーの調査========================
#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
if [ -f $plateau_width_output_file1 ]; then
    echo "Filename : ${plateau_width_output_file1}" 
    echo "File already exists. Please enter a new filename({filename}.csv):"
    read new_filename
    plateau_width_output_file1="${new_filename}.csv"
fi

for file in $(ls ../output/Jg04/*.csv); do
    input_file="$file" #入力ファイル名
    data=$(awk 'NR == 5 { print $3, $6 }' $input_file)
    echo $data >> $plateau_width_output_file1
done


#=============================5/9プラトーの調査========================
#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
if [ -f $plateau_width_output_file2 ]; then
    echo "Filename : ${plateau_width_output_file2}" 
    echo "File already exists. Please enter a new filename({filename}.csv):"
    read new_filename
    plateau_width_output_file2="${new_filename}.csv"
fi

for file in $(ls ../output/Jg04/*.csv); do
    input_file="$file" #入力ファイル名
    data=$(awk 'NR == 8 { print $3, $6 }' $input_file)
    echo $data >> $plateau_width_output_file2
done

#=============================7/9プラトーの調査========================
#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
if [ -f $plateau_width_output_file3 ]; then
    echo "Filename : ${plateau_width_output_file3}" 
    echo "File already exists. Please enter a new filename({filename}.csv):"
    read new_filename
    plateau_width_output_file3="${new_filename}.csv"
fi

for file in $(ls ../output/Jg04/*.csv); do
    input_file="$file" #入力ファイル名
    data=$(awk 'NR == 11 { print $3, $6 }' $input_file)
    echo $data >> $plateau_width_output_file3
done
