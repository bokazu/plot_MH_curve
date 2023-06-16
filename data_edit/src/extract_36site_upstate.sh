#!/bin/bash

input_dir="$1"
input_file=$(basename "${input_dir%.*}")
output_file="../output/${input_file}_upstate.csv"


for i in $(seq 1 13); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==9 { print }' "$input_dir" >> "$output_file"

for i in $(seq 15 23); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==18 { print }' "$input_dir" >> "$output_file"

for i in $(seq 24 32); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==27 { print }' "$input_dir" >> "$output_file"

for i in $(seq 33 34); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==1 { print }' "$input_dir" >> "$output_file"
awk 'NR==5 { print }' "$input_dir" >> "$output_file"


for i in $(seq 35 36); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==2 { print }' "$input_dir" >> "$output_file"
awk 'NR==6 { print }' "$input_dir" >> "$output_file"
