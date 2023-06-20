include("./src/kagome_spinrel.jl")
using .kagome_spinrel
using Plots
using Measures
using LaTeXStrings
anim = Animation()

for i=18:36
    input_szz_filename = "settings/36site/szz_" * string(i) * "_nn_list.csv"
    input_sz_filename = "settings/36site/sz_"* string(i) *"_upstate.csv"
    marker_magnification = 15
    line_magnification = 20
    output_szz_filename = "img/szz_rel_36site_"*string(i)*".png"
    output_szz_c_filename = "img/szz_c_rel_36site_"*string(i)*".png"
    plt_szz, plt_szz_c = kagome_spinrel.plot_36site_szz_NNcorellation(input_szz_filename,input_sz_filename ,marker_magnification,line_magnification,output_szz_filename,output_szz_c_filename)

    input_sxx_filename = "settings/36site/sxx_" * string(i) * "_nn_list.csv"
    input_sz_filename = "settings/36site/sz_"* string(i) *"_upstate.csv"
    marker_magnification = 15
    line_magnification = 20
    output_sxx_filename = "./img/sxx_rel_"*string(i)*".png"
    plt_sxx = kagome_spinrel.plot_36site_sxx_NNcorellation(input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_sxx_filename)
    
    MHdata_file="settings/MHdata.csv"
    MHdata_output="img/MHcurve.png"
    plt_MHcurve = kagome_spinrel.plot_36site_MH_curve_with_Highrigt(MHdata_file, i)

    M=[L"0",L"1/18",L"1/9",L"1/6",L"2/9",L"5/18",L"1/3 ",L"7/18",L"4/9",L"1/2",L"5/9",L"11/18",L"2/3",L"13/18",L"7/9 ",L"15/18",L"8/9 ",L"17/18",L"1  ",L"1"]
    l = @layout[a{0.01h}; b c d]
    title=plot(title="\n\n"*L"$\mathrm{Y}_3\mathrm{Cu}_9(\mathrm{OH})_{19}\mathrm{Cl}_8\ \ ,\ \ 36\mathrm{site}\ \ ,\ \  M/M_{sat}=$" * M[i-13],grid=false, titleposition =:left,showaxis=false,titleframe=:box)
    plt_summary = plot(title ,plt_MHcurve, plt_szz_c,plt_sxx, layout=l,size=(1640,640),margin=10mm, plot_titlevspan=0.01)
    frame(anim, plt_summary)
end

gif(anim, "36site_MHcurve_Ykapellasite.gif",fps=1)
