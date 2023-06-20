#!/bin/bash

szz_input_file_list=(../settings/szz_*)
sxx_input_file_list=(../settings/sxx_*)
sz_input_file_list=(../settings/sz_*)

rm ../output/*

for file in "${szz_input_file_list[@]}"; do
    echo "input file : $file"
    . extract_nn_36sites.sh "$file"
done

for file in "${sxx_input_file_list[@]}"; do
    echo "input file : $file"
    . extract_nn_36sites.sh "$file"
done

for file in "${sz_input_file_list[@]}"; do
    . extract_36sites_upstate.sh "$file"
done
