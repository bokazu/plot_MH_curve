include("./src/kagome_spinrel.jl")
using .kagome_spinrel
using Plots; gr()
using Measures
using LaTeXStrings

#格子点とbondの強調具合の指定
marker_magnification = 20
line_magnification = 20
#入力ファイル、出力ファイルの存在するディレクトリの指定
input_dir_name = "./settings/27site/real"
output_dir_name = "./img/27site/real"
#各ボンドのstrength
J_r = 0.87
J_g = 0.06
J_b = 1.0
# plt=kagome_spinrel.plot_27site_kagome_lattice(J_r,J_g,J_b)
# savefig(plt, "27site_lattice_of_Ykapllasite.png")

title_name = "\n\n" *L"$ \mathrm{Uniform} \ \ ,\ \  M/M_{sat}=$"
title_name = "\n\n"*L"$\mathrm{Y}_3\mathrm{Cu}_9(\mathrm{OH})_{19}\mathrm{Cl}_8\ \ ,\ \ 27\mathrm{site}\ \ ,\ \  , J_r = " * string(J_r) * "J_g = " * string(J_g) * "J_b = " * string(J_b) * L"$\ \ M/M_{sat}=$"
# 最低磁化から飽和磁化までの磁化曲線及び、各磁化でのスピン相関の様子をplotする
kagome_spinrel.summary_plot_27site_kagome(marker_magnification, line_magnification, input_dir_name, output_dir_name, J_r, J_g, J_b, title_name)