#!/bin/bash

merge_file1="eigen_val_plot_1.csv"
merge_file2="eigen_val_plot_2.csv"


#merge_eval_data.cppの実行
g++ merge_eval_data.cpp -o merge_eval_data
merge_eval_data 