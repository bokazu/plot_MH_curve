include("./src/kagome_spinrel.jl")
using .kagome_spinrel
using Plots             
using Measures
using LaTeXStrings

#格子点とbondの強調具合の指定
marker_magnification = 13
line_magnification = 13
# #入力ファイル、出力ファイルの存在するディレクトリの指定 注意)最後に'/'はつけない!
input_dir_name = "./settings/36site/Ykapellasite/real"
output_dir_name = "./img/Ykapellasite/36site/real"
# #各ボンドのstrength
J_r = 0.87
J_g = 0.06
J_b = 1.0

M_start = 18
M_end = 36
# title_name = "\n\n" *L"$ \mathrm{Uniform} \ \ ,\ \  M/M_{sat}=$"
title_name = "\n\n"*L"$\mathrm{Y}_3\mathrm{Cu}_9(\mathrm{OH})_{19}\mathrm{Cl}_8\ \ ,\ \ 36\mathrm{site}\ \ \ \  , J_r = $" * string(J_r) * L"$, J_g = $" * string(J_g) * L"$, J_b = $" * string(J_b) * L"$ \ \  , M/M_{sat}= $"
# #最低磁化から飽和磁化までの磁化曲線及び、各磁化でのスピン相関の様子をplotする
kagome_spinrel.summary_plot_36site_kagome(M_start,M_end,marker_magnification, line_magnification, input_dir_name, output_dir_name, J_r, J_g, J_b, title_name)


#===================-格子だけのplot===============================#
#格子点とbondの強調具合の指定
# marker_magnification = 20
# line_magnification = 20
# #入力ファイル、出力ファイルの存在するディレクトリの指定
# input_dir_name = "./settings/36site/real"
# output_dir_name = "./img/36site/real"
# #各ボンドのstrength
# J_r = 0.87
# J_g = 0.06
# J_b = 1.0

# plt=kagome_spinrel.plot_36site_kagome_lattice(J_r,J_g,J_b)
# savefig(plt, "36site_lattice_of_Ykapllasite.png")