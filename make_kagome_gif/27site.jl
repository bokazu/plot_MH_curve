include("./src/kagome_spinrel.jl")
using .kagome_spinrel
using Plots
using Measures
using LaTeXStrings
anim = Animation()

for i=14:27
    input_szz_filename = "settings/27site/uniform/szz_" * string(i) * "_nn_list.csv"
    input_sz_filename = "settings/27site/uniform/sz_"* string(i) *"_upstate.csv"
    marker_magnification = 10
    line_magnification = 30
    output_szz_filename = "img/szz_rel_"*string(i)*".png"
    output_szz_c_filename = "img/szz_c_rel_"*string(i)*".png"
    plt_szz, plt_szz_c = kagome_spinrel.plot_27site_szz_NNcorellation(input_szz_filename,input_sz_filename ,marker_magnification,line_magnification,output_szz_filename,output_szz_c_filename)

    input_sxx_filename = "settings/27site/uniform/sxx_" * string(i) * "_nn_list.csv"
    input_sz_filename = "settings/27site/uniform/sz_"* string(i) *"_upstate.csv"
    marker_magnification = 10
    line_magnification = 30
    output_sxx_filename = "./img/sxx_rel_"*string(i)*".png"
    plt_sxx = kagome_spinrel.plot_27site_sxx_NNcorellation(input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_sxx_filename)
    
    MHdata_file="settings/27site/uniform/MHdata.csv"
    MHdata_output="img/MHcurve.png"
    plt_MHcurve = kagome_spinrel.plot_27site_MH_curve_with_Highrigt(MHdata_file, i)

    M=[L"1/27",L"1/9",L"5/27",L"7/27",L"1/3 ",L"11/27",L"13/27",L"5/9 ",L"17/27",L"19/27",L"7/9 ",L"23/27 ",L"25/27",L"1  ",L"1  "]
    l = @layout[a{0.01h}; b c d]
    title=plot(title="\n\n" *L"$ \mathrm{Uniform} \ \ ,\ \  M/M_{sat}=$" * M[i-13],grid=false, titleposition =:left,showaxis=false,titleframe=:box)
    # title=plot(title="\n\n"*L"$\mathrm{Y}_3\mathrm{Cu}_9(\mathrm{OH})_{19}\mathrm{Cl}_8\ \ ,\ \ 27\mathrm{site}\ \ ,\ \  M/M_{sat}=$" * M[i-13],grid=false, titleposition =:left,showaxis=false,titleframe=:box)
    plt_summary = plot(title ,plt_MHcurve, plt_szz_c,plt_sxx, layout=l,size=(1640,640),margin=10mm, plot_titlevspan=0.01)
    frame(anim, plt_summary)
end

gif(anim, "27site_MHcurve_uniform.gif",fps=1)