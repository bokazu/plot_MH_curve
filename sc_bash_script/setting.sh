#!/bin/bash

# ===================並列化の設定===========================                                         
#SBATCH -p F1cpu                                                                                     
#SBATCH -N 1                                                                                         
#SBATCH -n 1                                                                                         
#SBATCH -c 128                                                                                       
#SBATCH --job-name="MHdata2(distorted)"                                                              
#SBATCH --mail-type=BEGIN                                                                            
#SBATCH --mail-type=END                                                                              
#SBATCH --mail-user=6222530@ed.tus.ac.jp 

ParamNum=20
JOB_NUM=$((2 * ParamNum)) #36site系の計算ではこちらを利用する
LATTICE="kagome"
SITENUM="36"
Lanczos_type="N" #lanczos法で固有値のみ求める => N /固有ベクトルも求める => V

#相互作用の大小関係は J_green < J_red < J_blue = 1
J_red=0.84
J_green=0.01
J_blue=1 #J_blueは1で固定
J_red_inc_steps=0.05 

var1="" #系の総サイト数
var2="" #部分系2つ+部分系をつなぐbondの本数(sys_numに対応)
var3="" #部分系1のサイト数(sys_site_Aに対応)
var4="" #部分系2のサイト数(sys_site_Bに対応)
var5="" #調べる部分空間について、磁化の最小値を設定する(min_up_spinに対応) 
var6="" #調べる部分空間について、磁化の最大値を設定する(max_up_spinに対応)
var7="" #Y-kapellasiteの模型におけるbondの相互作用J_red
var8="" #Y-kapellasiteの模型におけるbondの相互作用J_green
var9="" #Y-kapellasiteの模型におけるbondの相互作用J_blue
dir_jset0=""
dir_jset1=""
dir_jset2=""

#既存のsettingsディレクトリを削除する
rm -r settings_param*
#前回の計算結果を削除する
rm -r ./output/${LATTICE}/${SITENUM}site/*
#Jsetファイルを削除
rm -r ./sample_lists/${LATTICE}/${SITENUM}site/!(generate_jset|generate_jset.cpp|generate_jset.sh)


for p in $(seq 1 ${ParamNum});do
    #=======================Jsetファイルを用意する==================================
    dir_jset_output="./sample_lists/${LATTICE}/${SITENUM}site/param${p}"
    mkdir $dir_jset_output
    #bondの値は環境変数にしたので、以下のスクリプト内で上記の変数を読み込む
    ./sample_lists/${LATTICE}/${SITENUM}site/generate_jset "$J_red" "$J_green" "$J_blue" "$dir_jset_output"
    J_red=$(echo "$J_red + $J_red_inc_steps" | bc) #どのbondの値を変化させるかはここで決める
    #==================================================================================


    #=====outputディレクトリ内に各jobの結果を出力するためのディレクトリを作成する=======
    output_dir_base=$(basename "$dir_jset_output")
    
    dir_output_sxx_rel="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/sxx_rel"
    dir_output_szz_rel="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/szz_rel"
    dir_output_eval="./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}/eigen_val_"

    mkdir ./output/${LATTICE}/${SITENUM}/output_${output_dir_base}
    mkdir ./output/${LATTICE}/${SITENUM}/output_${output_dir_base}/sxx_rel
    mkdir ./output/${LATTICE}/${SITENUM}/output_${output_dir_base}/szz_rel
    mkdir ./output/${LATTICE}/${SITENUM}/output_${output_dir_base}/time
    #=====================================================================================
    

    #===============settingsディレクトリを用意しJsetファイルをコピーする===============
    for job in $(seq 1 ${JOB_NUM});do
        settings_dir="./settings_param${p}_job${job}"
        mkdir $settings_dir
        cp "$dir_jset_output"/* $settings_dir
    done
    #==================================================================================

    #system_info.txtはoutputディレクトリにコピーしておく
    cp $dir_jset_output/system_info.txt ./output/${LATTICE}/${SITENUM}site/output_${output_dir_base}

    for job in $(seq 1 ${JOB_NUM});do
        #system info.txtから系のサイト数などの情報を読み取るための処理
        for j in $(seq 1 9); do
        var="var$j"
        read -r $var < <(sed "${j}q;d" ./settings_param${p}_job${job}/system_info.txt)
        done

        #実行ファイルの作成
        PLOT_MH_EXE_FILE="main_param${p}_job${job}"
        export PLOT_MH_EXE_FILE
        cmake --build build

        export KMP_AFFINITY=scatter
        export OMP_SCHEDULE="dynamic,3"
        if (($job % 2 == 0)); then
            start_up_spin=$var5
            end_up_spin=$((var5 + 1))
            #[ToDo]jsetのディレクトリも渡す！
            sbatch ./F1cpu.sh "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_jset0" "$dir_jset1" "$dir_jset2" "$dir_output_eval" "$dir_output_time" "$dir_output_sxx_rel" "$dir_output_szz_rel" "$start_up_spin" "$end_up_spin" "$Lanczos_type"
        else
            start_up_spin=$((var5 + 2))
            end_up_spin=$var6
            #[ToDo]jsetのディレクトリも渡す！
            sbatch ./F1cpu.sh "$var2" "$var3" "$var4" "$var5" "$var6" "$var7" "$var8" "$var9" "$dir_jset0" "$dir_jset1" "$dir_jset2" "$dir_output_eval" "$dir_output_time" "$dir_output_sxx_rel" "$dir_output_szz_rel" "$start_up_spin" "$end_up_spin" "$Lanczos_type"
    done
done

scp k038112@ohtaka.issp.u-tokyo.ac.jp:output/kagome/27site/ .
