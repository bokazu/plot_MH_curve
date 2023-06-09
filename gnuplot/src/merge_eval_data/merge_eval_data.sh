#!/bin/bash

eval_data_file1=""
eval_data_file2=""
output_eval_data_file=""

echo "merge file 1 : $eval_data_file1"
echo "merge file 2 : $eval_data_file2"
echo "output file  : $output_eval_data_file" 

touch $output_eval_data_file
cat $eval_data_file1 > $output_eval_data_file
tail -n +1 $eval_data_file2 > $output_eval_data_file