include("./src/kagome_spinrel.jl")
using .kagome_spinrel
using Plots; gr()
using Measures
using LaTeXStrings
# using WebIO



# anim = Animation()
# for i=14:27
    # input_szz_filename = "settings/27site/distorted_1_0_1/szz_" * string(i) * "_nn_list.csv"
    # input_sz_filename = "settings/27site/distorted_1_0_1/sz_"* string(i) *"_upstate.csv"
    # marker_magnification = 20
    # line_magnification = 40
    # output_szz_filename = "img/27site/distorted_1_0_1/szz_rel_" * string(i) * ".png"
    # output_szz_c_filename = "img/27site/distorted_1_0_1/szz_c_rel_"*string(i)*".png"
    # plt_szz, plt_szz_c = kagome_spinrel.plot_27site_szz_NNcorellation(input_szz_filename,input_sz_filename ,marker_magnification,line_magnification,output_szz_filename,output_szz_c_filename)
    # input_sxx_filename = "settings/27site/distorted_1_0_1/sxx_" * string(i) * "_nn_list.csv"
    # input_sz_filename = "settings/27site/distorted_1_0_1/sz_"* string(i) *"_upstate.csv"
    # marker_magnification = 20
    # line_magnification = 40
    # output_sxx_filename = "./img/27site/distorted_1_0_1/sxx_rel_"*string(i)*".png"
    # plt_sxx = kagome_spinrel.plot_27site_sxx_NNcorellation(input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_sxx_filename)
    
    # MHdata_file="settings/27site/distorted_1_0_1/MHdata.csv"
    # MHdata_output="img/27site/distorted_1_0_1/MHcurve.png"
    # plt_MHcurve = kagome_spinrel.plot_27site_MH_curve_with_Highrigt(MHdata_file, i)

    # plt_lattice = kagome_spinrel.plot_27site_kagome_lattice(1, 0.0, 1.0)

    # M=[L"1/27",L"1/9",L"5/27",L"7/27",L"1/3 ",L"11/27",L"13/27",L"5/9 ",L"17/27",L"19/27",L"7/9 ",L"23/27 ",L"25/27",L"1  ",L"1  "]
    # l = @layout[a{0.01h}; b c; d e]
    # # title=plot(title="\n\n" *L"$ \mathrm{Uniform} \ \ ,\ \  M/M_{sat}=$" * M[i-13],grid=false, titleposition =:left,showaxis=false,titleframe=:box)
    # title=plot(title="\n\n"*L"$\mathrm{Y}_3\mathrm{Cu}_9(\mathrm{OH})_{19}\mathrm{Cl}_8\ \ ,\ \ 36\mathrm{site}\ \ ,\ \  M/M_{sat}=$" * M[i-13],grid=false, titleposition =:left,showaxis=false,titleframe=:box)
    # plt_summary = plot(title ,plt_MHcurve, plt_lattice,plt_szz_c,plt_sxx, layout=l,size=(1640,1640),margin=10mm, plot_titlevspan=0.01)
    # savefig("./img/27site/distorted_1_0_1/summary/summary_$i.pdf")
    # frame(anim, plt_summary)
# end

# gif(anim, "27site_MHcurve_distorted_1_0_1.gif",fps=1)

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