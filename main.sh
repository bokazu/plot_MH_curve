#!/bin/bash

target_dir="./sample_lists/kagome/27site/param6/*.txt"

var1="" #系の総サイト数
var2="" #部分系2つ+部分系をつなぐbondの本数(sys_numに対応)
var3="" #部分系1のサイト数(sys_site_Aに対応)
var4="" #部分系2のサイト数(sys_site_Bに対応)
var5="" #調べる部分空間について、磁化の最小値を設定する(min_up_spinに対応) 
var6="" #調べる部分空間について、磁化の最大値を設定する(max_up_spinに対応)

cp $target_dir ./settings
 
for i in $(seq 1 6); do
        var="var$i"
        read -r $var < <(sed "${i}q;d" ./settings/system_info.txt)
done


dir_output="./output/data_0.txt"

#===================コードを実行する==============================
cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
cmake --build build
./build/main_app "$var2" "$var3" "$var4" "$var5" "$var6" "$dir_output"

echo "Deleting files in ./settings/"
rm ./settings/*.txt