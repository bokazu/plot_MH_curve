#!/bin/bash
ParamNum=1
LATTICE="kagome"
SITENUM="27"
Lanczos_type="N" #lanczos法で固有値のみ求める => N /固有ベクトルも求める => V
terminal_output_file="./terminal_output_${LATTICE}_${SITENUM}.txt"
touch $terminal_output_file

#相互作用の大小関係は J_green < J_red < J_blue = 1
J_red=1.0
J_green=1.0
J_blue=1.0 #J_blueは1で固定
J_red_inc_steps=0.2

#実行ファイル名の指定
PLOT_MH_EXE_FILE="main"

#build
export PLOT_MH_EXE_FILE

cmake -S . -DCMAKE_CXX_COMPILER=icpx -B build
cmake --build build | tee $terminal_output_file

var1="" #系の総サイト数
var2="" #部分系2つ+部分系をつなぐbondの本数(sys_numに対応)
var3="" #部分系1のサイト数(sys_site_Aに対応)
var4="" #部分系2のサイト数(sys_site_Bに対応)
var5="" #調べる部分空間について、磁化の最小値を設定する(min_up_spinに対応) 
var6="" #調べる部分空間について、磁化の最大値を設定する(max_up_spinに対応)
var7="" #Y-kapellasiteの模型におけるbondの相互作用J_red
var8="" #Y-kapellasiteの模型におけるbondの相互作用J_green
var9="" #Y-kapellasiteの模型におけるbondの相互作用J_blue


#前回の計算結果を削除
rm -r ./output/${LATTICE}/${SITENUM}site/!(real|hexagon|trimer|jred)  | tee -a $terminal_output_file
#Jsetファイルを削除
rm -r ./sample_lists/${LATTICE}/${SITENUM}site/!(generate_jset|generate_jset.cpp|generate_jset.sh)  | tee -a $terminal_output_file

for p in $(seq 1 ${ParamNum});do
    #=======================Jsetファイルを用意する==================================
    dir_jset_output="./sample_lists/${LATTICE}/${SITENUM}site/param${p}"
    mkdir $dir_jset_output | tee -a $terminal_output_file
    #bondの値は環境変数にしたので、以下のスクリプト内で上記の変数を読み込む
    ./sample_lists/${LATTICE}/${SITENUM}site/generate_jset "$J_red" "$J_green" "$J_blue" "$dir_jset_output" | tee -a $terminal_output_file
    J_red=$(echo "$J_red + $J_red_inc_steps" | bc) #どのbondの値を変化させるかはここで決める
    #==================================================================================

    #===================Jsetファイルをsettingsにコピーする============================
    cp "$dir_jset_output"/* ./settings/ | tee -a $terminal_output_file
    dir_jset0="./settings/jset0.txt"
    dir_jset1="./settings/jset1.txt"
    dir_jset2="./settings/jset2.txt"
    #==================================================================================


    #=====outputディレクトリ内に各jobの結果を出力するためのディレクトリを作成する=======
    output_dir_base=$(basename "$dir_jset_output")
    
    dir_output_sxx_rel="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/sxx_rel"
    dir_output_sz_rel="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/sz_rel"
    dir_output_szz_rel="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/szz_rel"
    dir_output_time="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/time"
    file_output_eval="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/eigen_val.csv"
    file_output_time="$dir_output_time/time_info_"

    mkdir ./output/${LATTICE}/${SITENUM}site/output_${output_dir_base} | tee -a $terminal_output_file
    mkdir $dir_output_sxx_rel | tee -a $terminal_output_file
    mkdir $dir_output_szz_rel | tee -a $terminal_output_file
    mkdir $dir_output_sz_rel | tee -a $terminal_output_file
    mkdir $dir_output_time | tee -a $terminal_output_file
    #=====================================================================================

    #system_info.txtはoutputディレクトリにコピーしておく
    cp $dir_jset_output/system_info.txt ./output/${LATTICE}/${SITENUM}site/output_${output_dir_base} | tee -a $terminal_output_file

    #system info.txtから系のサイト数などの情報を読み取るための処理
    for j in $(seq 1 9); do
        var="var$j"
        read -r $var < <(sed "${j}q;d" ./settings/system_info.txt) #ここにteeコマンドを入れてはいけない
        echo $var
    done

    export OMP_NUM_THREADS=1 #変更箇所2
    # export KMP_AFFINITY=scatter
    # export OMP_SCHEDULE="dynamic,3"

    start_up_spin=$var5
    end_up_spin=$var5 #変更箇所1
    ./build/${PLOT_MH_EXE_FILE} "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_jset0" "$dir_jset1" "$dir_jset2" "$file_output_eval" "$file_output_time" "$dir_output_sxx_rel" "$dir_output_sz_rel" "$dir_output_szz_rel" "$start_up_spin" "$end_up_spin" "$Lanczos_type" | tee -a $terminal_output_file
    
    #===================jsetファイルをsettings_dirから削除する==============
    echo "Deleting files in ./settings" | tee -a $terminal_output_file
    rm ./settings/* | tee -a $terminal_output_file
done
