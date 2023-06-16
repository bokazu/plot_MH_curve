#!/bin/bash

input_dir="$1"
input_file=$(basename "${input_dir%.*}")
output_file="../output/${input_file}_upstate.csv"


for i in $(seq 1 9); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==4 { print }' "$input_dir" >> "$output_file"

for i in $(seq 10 12); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

for i in $(seq 13 18); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==13 { print }' "$input_dir" >> "$output_file"

for i in $(seq 19 27); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done

awk 'NR==22 { print }' "$input_dir" >> "$output_file"

for i in $(seq 1 3); do
    awk -v line=$i 'NR==line { print }' "$input_dir" >> "$output_file"
done