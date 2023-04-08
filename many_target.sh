#!/bin/bash

#カゴメ格子27site系の情報
site="27"
min_up_spin="14"
max_up_spin="27"
start_param=0

#計算を行いたいjsetファイルを格納しているディレクトリのリストを取得する
from_dir="./sample_lists/kagome/27site"
to_dir="./settings"

dir_list=("$from_dir"/*/)

cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build


for dir in "${dir_list[@]}"; do
    dir=${dir%*/}
    echo "Copying files from $dir to $to_dir"
    cp -v "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

    #===================コードを実行する==============================
    #----------------------------変更箇所--------------------------#
    sys_num="14"
    sys_site_A="18"
    sys_site_B="9"
    dir_output="./output/data_${start_param}.txt"

    cmake --build build
    ./build/main_app "$sys_num" "$sys_site_A" "$sys_site_B" "$min_up_spin" "$max_up_spin" "$dir_output"

    start_param=$((start_param + 1))
    #-----------------jsetファイルをto_dirから削除する-----------------
    echo "Deleting files in $to_dir"
    rm -v $to_dir/*
    
done