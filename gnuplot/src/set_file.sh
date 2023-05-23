#!/bin/bash


#===================相互作用とplateau幅の関係を調べるための処理===
J_num="Jg07"

plateau_width_output_file1="../../output/plateau_width/${J_num}/plateau_width_dif_${J_num}_3of9.csv" #[コメント]projectタイトルのような変数を用意して、それを末尾にくっつけると良いかも?
plateau_width_output_file2="../../output/plateau_width/${J_num}/plateau_width_dif_${J_num}_5of9.csv"
plateau_width_output_file3="../../output/plateau_width/${J_num}/plateau_width_dif_${J_num}_7of9.csv"

plateau_width_output_file4="../../output/plateau_width/${J_num}/plateau_width_dif_${J_num}_1of9.csv"
#=============================入力データが存在するかチェック========================
#同一ファイル名のファイルが存在しないかを確認し、存在する場合にはファイル名を変更する
#1/3プラトー用入力ファイルのチェック
# if [ -f $plateau_width_output_file1 ]; then
#     echo "Filename : ${plateau_width_output_file1}" 
#     echo "File already exists. Please enter a new filename(../../output/plateau_width/${J_num}/{filename}.csv}):"
#     read new_filename
#     plateau_width_output_file1="../../output/plateau_width/${J_num}/${new_filename}.csv"
# fi

# #5/9プラトー用入力ファイルのチェック
# if [ -f $plateau_width_output_file2 ]; then
#     echo "Filename : ${plateau_width_output_file2}" 
#     echo "File already exists. Please enter a new filename(../../output/plateau_width/${J_num}/{filename}.csv):"
#     read new_filename
#     plateau_width_output_file2="../../output/plateau_width/${J_num}/${new_filename}.csv"
# fi

# #7/9プラトー用入力ファイルのチェック
# if [ -f $plateau_width_output_file3 ]; then
#     echo "Filename : ${plateau_width_output_file3}" 
#     echo "File already exists. Please enter a new filename(../../output/{filename}.csv):"
#     read new_filename
#     plateau_width_output_file3="../../output/${new_filename}.csv"
# fi

#入力ファイルについての処理
for file in $(ls ../../output/${J_num}/*.csv); do
    input_file="$file" #入力ファイル名
    
    #1/9プラトーのデータ取得
    data_1of9=$(awk 'NR == 2 { print $3, $6 }' $input_file)
    echo $data_1of9 >> $plateau_width_output_file4

    #1/3プラトーのデータ取得
    data_1of3=$(awk 'NR == 5 { print $3, $6 }' $input_file)
    echo $data_1of3 >> $plateau_width_output_file1

    #5/9プラトーのデータ取得
    data_5of9=$(awk 'NR == 8 { print $3, $6 }' $input_file)
    echo $data_5of9 >> $plateau_width_output_file2

    #7/9プラトーのデータ取得
    data_7of9=$(awk 'NR == 11 { print $3, $6 }' $input_file)
    echo $data_7of9 >> $plateau_width_output_file3
done
