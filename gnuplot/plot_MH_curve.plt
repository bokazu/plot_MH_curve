set key left top
set xrange[0:3.5]
plot "../output/MHdata_0.csv" using 1:2 with steps title "param0"
replot "../output/MHdata_5.csv" using 1:2 with steps title "param5"
replot "../output/MHdata_6.csv" using 1:2 with steps title "param6"
replot "../output/MHdata_7.csv" using 1:2 with steps title "param7"

set title "kagome 27 site various J_green"

set term png
set output "kagome_27site_various_J_green.png";replot
set term qt