include("./src/kagome_spinrel.jl")
using .kagome_spinrel
using Plots
using Measures
anim = Animation()

for i=14:27
    input_szz_filename = "settings/szz_" * string(i) * "_nn_list.csv"
    input_sz_filename = "settings/sz_"* string(i) *"_upstate.csv"
    marker_magnification = 50
    line_magnification = 20
    output_szz_filename = "img/szz_rel_"*string(i)*".png"
    output_szz_c_filename = "img/szz_c_rel_"*string(i)*".png"
    plt_szz, plt_szz_c = kagome_spinrel.plot_27site_szz_NNcorellation(input_szz_filename,input_sz_filename ,marker_magnification,line_magnification,output_szz_filename,output_szz_c_filename)

    input_sxx_filename = "settings/sxx_" * string(i) * "_nn_list.csv"
    input_sz_filename = "settings/sz_"* string(i) *"_upstate.csv"
    marker_magnification = 50
    line_magnification = 20
    output_sxx_filename = "./img/sxx_rel_"*string(i)*".png"
    plt_sxx = kagome_spinrel.plot_27site_sxx_NNcorellation(input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_sxx_filename)
    
    MHdata_file="settings/MHdata.csv"
    MHdata_output="img/MHcurve.png"
    plt_MHcurve = kagome_spinrel.plot_27site_MH_curve_with_Highrigt(MHdata_file, i)

    M=["1/27","1/9","5/27","7/27","1/3","11/27","13/27","5/9","17/27","19/27","7/9","23/27","25/27","1","1"]
    l = @layout[a{0.01h}; b c d]
    title=plot(title="M/M_sat=" * M[i-13],grid=false, showaxis=false)
    plt_summary = plot(title ,plt_MHcurve, plt_szz_c,plt_sxx, layout=l,size=(1640,560),margin=5mm)

    frame(anim, plt_summary)
end

gif(anim, "MHcurve_modify.gif",fps=1)
