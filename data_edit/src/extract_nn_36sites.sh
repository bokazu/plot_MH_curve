#!/bin/bash

input_dir="$1"
input_file=$(basename "${input_dir%.*}")
output_file="../output/${input_file}_nn_list.csv"

awk 'NR==3 { print }' "$input_dir" >> "$output_file"
awk 'NR==4 { print }' "$input_dir" >> "$output_file"
awk 'NR==16 { print }' "$input_dir" >> "$output_file"
awk 'NR==20 { print }' "$input_dir" >> "$output_file"
awk 'NR==29 { print }' "$input_dir" >> "$output_file"
awk 'NR==32 { print }' "$input_dir" >> "$output_file"
awk 'NR==41 { print }' "$input_dir" >> "$output_file"
awk 'NR==44 { print }' "$input_dir" >> "$output_file"
awk 'NR==54 { print }' "$input_dir" >> "$output_file"
awk 'NR==65 { print }' "$input_dir" >> "$output_file"
awk 'NR==66 { print }' "$input_dir" >> "$output_file"
awk 'NR==75 { print }' "$input_dir" >> "$output_file"
awk 'NR==76 { print }' "$input_dir" >> "$output_file"
awk 'NR==84 { print }' "$input_dir" >> "$output_file"
awk 'NR==79 { print }' "$input_dir" >> "$output_file"
awk 'NR==86 { print }' "$input_dir" >> "$output_file"
awk 'NR==535 { print }' "$input_dir" >> "$output_file"
awk 'NR==92 { print }' "$input_dir" >> "$output_file"
awk 'NR==557 { print }' "$input_dir" >> "$output_file"
awk 'NR==97 { print }' "$input_dir" >> "$output_file"
awk 'NR==580 { print }' "$input_dir" >> "$output_file"
awk 'NR==101 { print }' "$input_dir" >> "$output_file"
awk 'NR==602 { print }' "$input_dir" >> "$output_file"
awk 'NR==104 { print }' "$input_dir" >> "$output_file"
awk 'NR==625 { print }' "$input_dir" >> "$output_file"
awk 'NR==90 { print }' "$input_dir" >> "$output_file"
awk 'NR==647 { print }' "$input_dir" >> "$output_file"
awk 'NR==109 { print }' "$input_dir" >> "$output_file"
awk 'NR==110 { print }' "$input_dir" >> "$output_file"
awk 'NR==132 { print }' "$input_dir" >> "$output_file"
awk 'NR==133 { print }' "$input_dir" >> "$output_file"
awk 'NR==154 { print }' "$input_dir" >> "$output_file"
awk 'NR==155 { print }' "$input_dir" >> "$output_file"
awk 'NR==170 { print }' "$input_dir" >> "$output_file"
awk 'NR==189 { print }' "$input_dir" >> "$output_file"
awk 'NR==193 { print }' "$input_dir" >> "$output_file"
awk 'NR==207 { print }' "$input_dir" >> "$output_file"
awk 'NR==210 { print }' "$input_dir" >> "$output_file"
awk 'NR==224 { print }' "$input_dir" >> "$output_file"
awk 'NR==227 { print }' "$input_dir" >> "$output_file"
awk 'NR==240 { print }' "$input_dir" >> "$output_file"
awk 'NR==242 { print }' "$input_dir" >> "$output_file"
awk 'NR==174 { print }' "$input_dir" >> "$output_file"
awk 'NR==257 { print }' "$input_dir" >> "$output_file"
awk 'NR==177 { print }' "$input_dir" >> "$output_file"
awk 'NR==272 { print }' "$input_dir" >> "$output_file"
awk 'NR==273 { print }' "$input_dir" >> "$output_file"
awk 'NR==286 { print }' "$input_dir" >> "$output_file"
awk 'NR==287 { print }' "$input_dir" >> "$output_file"
awk 'NR==299 { print }' "$input_dir" >> "$output_file"
awk 'NR==294 { print }' "$input_dir" >> "$output_file"
awk 'NR==305 { print }' "$input_dir" >> "$output_file"
awk 'NR==310 { print }' "$input_dir" >> "$output_file"
awk 'NR==315 { print }' "$input_dir" >> "$output_file"
awk 'NR==319 { print }' "$input_dir" >> "$output_file"
awk 'NR==324 { print }' "$input_dir" >> "$output_file"
awk 'NR==328 { print }' "$input_dir" >> "$output_file"
awk 'NR==332 { print }' "$input_dir" >> "$output_file"
awk 'NR==335 { print }' "$input_dir" >> "$output_file"
awk 'NR==339 { print }' "$input_dir" >> "$output_file"
awk 'NR==375 { print }' "$input_dir" >> "$output_file"
awk 'NR==309 { print }' "$input_dir" >> "$output_file"
awk 'NR==376 { print }' "$input_dir" >> "$output_file"
awk 'NR==465 { print }' "$input_dir" >> "$output_file"
awk 'NR==351 { print }' "$input_dir" >> "$output_file"
awk 'NR==355 { print }' "$input_dir" >> "$output_file"
awk 'NR==400 { print }' "$input_dir" >> "$output_file"
awk 'NR==467 { print }' "$input_dir" >> "$output_file"
awk 'NR==357 { print }' "$input_dir" >> "$output_file"
awk 'NR==489 { print }' "$input_dir" >> "$output_file"
awk 'NR==402 { print }' "$input_dir" >> "$output_file"
awk 'NR==490 { print }' "$input_dir" >> "$output_file"