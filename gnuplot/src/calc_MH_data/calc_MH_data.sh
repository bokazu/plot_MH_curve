#!/bin/bash

eval_data_dir="../../../output/kagome/36site/output_param1"
# eval_data_dir="../../../output/kagome/36site/distorted/param3"
eval_data_file="${eval_data_dir}/eigen_val.csv"
MH_data_file="${eval_data_dir}/MHdata.csv"

# 格子の情報
var7="" #Y-kapellasiteの模型におけるbondの相互作用J_red
var8="" #Y-kapellasiteの模型におけるbondの相互作用J_green
var9="" #Y-kapellasiteの模型におけるbondの相互作用J_blue

#outputディレクトリにコピーしてきたsystem_info.txtから格子の情報を読み取る
for i in $(seq 7 9); do
    var="var$i"
    read -r $var < <(sed "${i}q;d" ${eval_data_dir}/system_info.txt)
done

#merge_eval_data.cppの実行
g++ calc_MH_data.cpp -o calc_MH_data 
./calc_MH_data "$var7" "$var8" "$var9" "$eval_data_file" "$MH_data_file"
