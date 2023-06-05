#!/bin/bash

# [Memo]これらの情報を各paramディレクトリ内のsystem_info.txtから読みこむためのコードを追加する
# カゴメ格子の情報
var1="" #系の総サイト数
var2="" #部分系2つ+部分系をつなぐbondの本数(sys_numに対応)
var3="" #部分系1のサイト数(sys_site_Aに対応)
var4="" #部分系2のサイト数(sys_site_Bに対応)
var5="" #調べる部分空間について、磁化の最小値を設定する(min_up_spinに対応) 
var6="" #調べる部分空間について、磁化の最大値を設定する(max_up_spinに対応)
var7="" #Y-kapellasiteの模型におけるbondの相互作用J_red
var8="" #Y-kapellasiteの模型におけるbondの相互作用J_green
var9="" #Y-kapellasiteの模型におけるbondの相互作用J_blue

LATTICE="kagome"
SITENUM="27"
# 計算を行いたいjsetファイルを格納しているディレクトリのリストを取得する
from_dir="./sample_lists/${LATTICE}/${SITENUM}site"
to_dir="./settings"
dir_list=("$from_dir"/*/)


cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
cmake --build build

export KMP_AFFINITY=scatter
export OMP_SCHEDULE="dynamic,3"

for dir in "${dir_list[@]}"; do
        dir=${dir%*/}
        echo "Copying files from $dir to $to_dir"
        cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

        #system_info.txtから系のサイト数などの情報を読み取るための処理
        for j in $(seq 1 9); do
            var="var$j"
            read -r $var < <(sed "${j}q;d" ./settings/system_info.txt)
        done

        param_dir=$(basename "$dir")
        mkdir ./output/${LATTICE}/${SITENUM}site/$param_dir
        mkdir ./output/${LATTICE}/${SITENUM}site/$param_dir/sxx_rel
        mkdir ./output/${LATTICE}/${SITENUM}site/$param_dir/szz_rel
        mkdir ./output/${LATTICE}/${SITENUM}site/$param_dir/time
        
        dir_output_sxx_rel="./output/${LATTICE}/${SITENUM}site/$param_dir/sxx_rel/spin_sxx_rel"
        dir_output_szz_rel="./output/${LATTICE}/${SITENUM}site/$param_dir/szz_rel/spin_szz_rel"
	dir_output_eval="./output/${LATTICE}/${SITENUM}site/$param_dir/eigen_val.csv"
	dir_output_time="./output/${LATTICE}/${SITENUM}site/$param_dir/time/time_info2_"

        #===================コードを実行する==============================           
        ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output_eval" "$dir_output_time" "$dir_output_sxx_rel" "$dir_output_szz_rel"

        #===================jsetファイルをto_dirから削除する==============
        echo "Deleting files in $to_dir"
        rm $to_dir/*
 done
