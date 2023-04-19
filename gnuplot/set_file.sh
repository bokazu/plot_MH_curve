#!/bin/bash


#===================相互作用とplateau幅の関係を調べるための処理===
plateau_width_output_file1="../output/plateau_width_dif_Jr_Jg_7of9.csv" #[コメント]projectタイトルのような変数を用意して、それを末尾にくっつけると良いかも?

#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
if [ -f $plateau_width_output_file1 ]; then
    echo "Filename : ${plateau_width_output_file1}" 
    echo "File already exists. Please enter a new filename({filename}.csv):"
    read new_filename
    plateau_width_output_file1="${new_filename}.csv"
fi


for file in $(ls ../output/*.csv); do
    input_file="$file" #入力ファイル名
    data=$(awk 'NR == 11 { print $3, $6 }' $input_file)
    echo $data >> $plateau_width_output_file1
done