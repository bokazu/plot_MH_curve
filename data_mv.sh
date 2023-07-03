#!/bin/bash

param_num=10
cp ./output/kagome/36site/distorted/jred/output_param${param_num}/sz_rel/sz_*.csv ./data_edit/settings/
cp ./output/kagome/36site/distorted/jred/output_param${param_num}/szz_rel/szz_*.csv ./data_edit/settings/
cp ./output/kagome/36site/distorted/jred/output_param${param_num}/sxx_rel/sxx_*.csv ./data_edit/settings/

cd ./data_edit/src
source data_edit_36site.sh
cd -
rm ./data_edit/settings/*.csv

mv ./data_edit/output/sz_*.csv ./make_kagome_gif/settings/36site/Ykapellasite/jred/param${param_num}/sz/
mv ./data_edit/output/szz_*.csv ./make_kagome_gif/settings/36site/Ykapellasite/jred/param${param_num}/szz/
mv ./data_edit/output/sxx_*.csv ./make_kagome_gif/settings/36site/Ykapellasite/jred/param${param_num}/sxx/
cp ./output/kagome/36site/distorted/jred/output_param${param_num}/system_info.txt ./make_kagome_gif/settings/36site/Ykapellasite/jred/param${param_num}/

