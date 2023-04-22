#!/bin/bash

# [Memo]これらの情報を各paramディレクトリ内のsystem_info.txtから読みこむためのコードを追加する
# カゴメ格子27site系の情報
var1="" #系の総サイト数
var2="" #部分系2つ+部分系をつなぐbondの本数(sys_numに対応)
var3="" #部分系1のサイト数(sys_site_Aに対応)
var4="" #部分系2のサイト数(sys_site_Bに対応)
var5="" #調べる部分空間について、磁化の最小値を設定する(min_up_spinに対応) 
var6="" #調べる部分空間について、磁化の最大値を設定する(max_up_spinに対応)
var7="" #Y-kapellasiteの模型におけるbondの相互作用J_red
var8="" #Y-kapellasiteの模型におけるbondの相互作用J_green
var9="" #Y-kapellasiteの模型におけるbondの相互作用J_blue

===================jsetファイルを用意する===========================



# 計算を行いたいjsetファイルを格納しているディレクトリのリストを取得する
from_dir="./sample_lists/kagome/27site"
to_dir="./settings"

dir_list=("$from_dir"/*/)

#./outputに入っているファイルを削除する
# rm -v ./output/*.txt ./output/*.csv


# for dir in "${dir_list[@]}"; do
#     dir=${dir%*/}
#     echo "Copying files from $dir to $to_dir"
#     cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

#     #system_info.txtから系のサイト数などの情報を読み取るための処理
#     for i in $(seq 1 9); do
#         var="var$i"
#         read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
#     done

#     param_dir=$(basename "$dir")
#     dir_output="./output/Jg03/MHdata_$param_dir.csv"

#     #===================コードを実行する==============================
#     cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
#     cmake --build build
#     ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

#     #===================jsetファイルをto_dirから削除する==============
#     echo "Deleting files in $to_dir"
#     rm $to_dir/*
# done
J_green=0.0
J_green_inc_steps=0.1

cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
dir_numbers=(00 01 02 03 04 05 06 07 08 09 10)
for dir_number in "${dir_numbers[@]}"; do
    cd sample_lists/kagome/27site/
    rm -r param*
    cd ../../

    . main_generate_jset.sh
    cd ..

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
        dir_output="./output/Jg${dir_number}/MHdata_$param_dir.csv"
    

        #===================コードを実行する==============================        
        cmake --build build
        ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

        #===================jsetファイルをto_dirから削除する==============
        echo "Deleting files in $to_dir"
        rm $to_dir/*
    done

    J_green=$(echo "$J_green + $J_green_inc_steps" | bc) #J_greenの値を更新する
done

# cd sample_lists/kagome/27site/
# rm -r param*
# cd ../../
# J_green=0.5
# . main_generate_jset.sh
# cd ..

# for dir in "${dir_list[@]}"; do
#     dir=${dir%*/}
#     echo "Copying files from $dir to $to_dir"
#     cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

#     #system_info.txtから系のサイト数などの情報を読み取るための処理
#     for i in $(seq 1 9); do
#         var="var$i"
#         read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
#     done

#     param_dir=$(basename "$dir")
#     dir_output="./output/Jg05/MHdata_$param_dir.csv"

#     #===================コードを実行する==============================
#     cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
#     cmake --build build
#     ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

#     #===================jsetファイルをto_dirから削除する==============
#     echo "Deleting files in $to_dir"
#     rm $to_dir/*
# done


# cd sample_lists/kagome/27site/
# rm -r param*
# cd ../../
# J_green=0.6
# . main_generate_jset.sh
# cd ..

# for dir in "${dir_list[@]}"; do
#     dir=${dir%*/}
#     echo "Copying files from $dir to $to_dir"
#     cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

#     #system_info.txtから系のサイト数などの情報を読み取るための処理
#     for i in $(seq 1 9); do
#         var="var$i"
#         read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
#     done

#     param_dir=$(basename "$dir")
#     dir_output="./output/Jg06/MHdata_$param_dir.csv"

#     #===================コードを実行する==============================
#     cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
#     cmake --build build
#     ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

#     #===================jsetファイルをto_dirから削除する==============
#     echo "Deleting files in $to_dir"
#     rm $to_dir/*
# done


# cd sample_lists/kagome/27site/
# rm -r param*
# cd ../../
# J_green=0.7
# . main_generate_jset.sh
# cd ..

# for dir in "${dir_list[@]}"; do
#     dir=${dir%*/}
#     echo "Copying files from $dir to $to_dir"
#     cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

#     #system_info.txtから系のサイト数などの情報を読み取るための処理
#     for i in $(seq 1 9); do
#         var="var$i"
#         read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
#     done

#     param_dir=$(basename "$dir")
#     dir_output="./output/Jg07/MHdata_$param_dir.csv"

#     #===================コードを実行する==============================
#     cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
#     cmake --build build
#     ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

#     #===================jsetファイルをto_dirから削除する==============
#     echo "Deleting files in $to_dir"
#     rm $to_dir/*
# done


# cd sample_lists/kagome/27site/
# rm -r param*
# cd ../../
# J_green=0.8
# . main_generate_jset.sh
# cd ..

# for dir in "${dir_list[@]}"; do
#     dir=${dir%*/}
#     echo "Copying files from $dir to $to_dir"
#     cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

#     #system_info.txtから系のサイト数などの情報を読み取るための処理
#     for i in $(seq 1 9); do
#         var="var$i"
#         read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
#     done

#     param_dir=$(basename "$dir")
#     dir_output="./output/Jg08/MHdata_$param_dir.csv"

#     #===================コードを実行する==============================
#     cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
#     cmake --build build
#     ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

#     #===================jsetファイルをto_dirから削除する==============
#     echo "Deleting files in $to_dir"
#     rm $to_dir/*
# done

# cd sample_lists/kagome/27site/
# rm -r param*
# cd ../../
# J_green=0.9
# . main_generate_jset.sh
# cd ..

# for dir in "${dir_list[@]}"; do
#     dir=${dir%*/}
#     echo "Copying files from $dir to $to_dir"
#     cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

#     #system_info.txtから系のサイト数などの情報を読み取るための処理
#     for i in $(seq 1 9); do
#         var="var$i"
#         read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
#     done

#     param_dir=$(basename "$dir")
#     dir_output="./output/Jg09/MHdata_$param_dir.csv"

#     #===================コードを実行する==============================
#     cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
#     cmake --build build
#     ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

#     #===================jsetファイルをto_dirから削除する==============
#     echo "Deleting files in $to_dir"
#     rm $to_dir/*
# done


# cd sample_lists/kagome/27site/
# rm -r param*
# cd ../../
# J_green=1.0
# . main_generate_jset.sh
# cd ..

# for dir in "${dir_list[@]}"; do
#     dir=${dir%*/}
#     echo "Copying files from $dir to $to_dir"
#     cp "$dir"/* $to_dir/ #jsetファイルを実行用のフォルダにコピーする

#     #system_info.txtから系のサイト数などの情報を読み取るための処理
#     for i in $(seq 1 9); do
#         var="var$i"
#         read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
#     done

#     param_dir=$(basename "$dir")
#     dir_output="./output/Jg10/MHdata_$param_dir.csv"

#     #===================コードを実行する==============================
#     cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
#     cmake --build build
#     ./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_output"

#     #===================jsetファイルをto_dirから削除する==============
#     echo "Deleting files in $to_dir"
#     rm $to_dir/*
# done
