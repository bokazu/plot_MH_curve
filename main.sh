#!/bin/bash


#========================================変更箇所=====================================#
site="27"
sys_num="14"
sys_site_A="18"
sys_site_B="9"
min_up_spin="14"
max_up_spin="27"
dir_output="output/data.txt"
parameter_num="0"

#sample_lists中のjsetテキストfileのコピー先。cppファイルはこのディレクトリ内のfileを使用する
run_dir="./settings"  #ここは変更しなくて良い

#bondの情報が書かれたテキストファイルをsettingsディレクトリにコピーしてくる
sample_dir="./sample_lists/kagome/${site}site/param${parameter_num}"
echo $sample_dir
#====================================================================================#

#上記ディレクトリ内にあるファイル名を代入する
sample_files=$(ls $sample_dir/jset*.txt)
echo $sample_files

#テキストファイルをコピーする
for f in $sample_files; do
    cp "$f" $run_dir
done

cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
cmake --build build

./build/main_app "$sys_num" "$sys_site_A" "$sys_site_B" "$min_up_spin" "$max_up_spin" "$dir_output"