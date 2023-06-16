

module kagome_spinrel
using Plots
using Measures
using LaTeXStrings

    function plot_27site_MH_curve_with_Highrigt(MHdata_file, num_of_M)
        h=zeros(15)
        M=zeros(15)
        #MHcurveのデータの読み込み
        open(MHdata_file) do MH_file
           MH_index=1
           for line in eachline(MH_file)
              s = split(line," ")
              h[MH_index]=parse(Float64,s[1])    
              M[MH_index]=parse(Float64,s[2])
              MH_index+=1
            end
        end

        h[15] = 100
        M[15] = 1.0
        plt = plot(h,M,line=:steppre,xlabel=L"h", ylabel=L"M/M_{sat}",xlims=(0,3.5),ylims=(0,1),title="Magnetization curve",titleposition=:center,legend=false, titlegap=2)
        if num_of_M < 26
            plot!(plt, [h[num_of_M-13],h[num_of_M-12]],[M[num_of_M-12],M[num_of_M-12]], lw=3,linecolor=:red, titlegap=2)
        elseif num_of_M == 26
            plot!(plt, [h[num_of_M-13],3.5],[M[num_of_M-12],M[num_of_M-12]],lw=3, linecolor=:red, titlegap=2)
        end
        return plt
    end


#<S_i^zS_j^z>,<S_i^zS_j^z>_Cについて最近接相互作用についてのみスピン相関の大きさを線の太さで符号を(色)でplotする。また、格子点は<S_i^z>を計算し、直径で大きさを、色で符号を表しplotする
#- input_szz_filename   : <S_i^zS_j^z>,<S_i^zS_j^z>_Cについて最近接サイトのデータのみが記載されたファイル名を指定する
#- input_sz_filename    : <S_i^z>のデータのみが記載されたファイル名を指定する
#- marker_magnification : 格子点の半径はデフォルトでは<S_i^z>の絶対値となるが、その値が小さすぎる場合に倍率をこの変数で指定する
#- line_magnification   : bond間の線の太さはデフォルトでは<S_i^zS_j^z>,<S_i^zS_j^z>_Cの絶対値となるが、その値が小さすぎる場合に倍率をこの変数で指定する
    function plot_27site_szz_NNcorellation(input_szz_filename,input_sz_filename ,marker_magnification,line_magnification,output_szz_filename,output_szz_c_filename)
        #1行目の格子点のx,y座標
        Y1 = 0
        x1=[0.0,3.0,6.0]
        y1=[0.0,0.0,0.0]

        #2行目の格子点のx,y座標
        Y2 = Y1-3*sqrt(3)/4
        x2=[-9/4,-3/4,3/4,9/4,15/4,21/4,27/4]
        y2=[Y2,    Y2,  Y2,Y2,  Y2,  Y2,  Y2]

        #3行目の格子点のx,y座標
        Y3 = Y2-3*sqrt(3)/4
        x3=[-3/2, 3/2,9/2]
        y3=[Y3,Y3,Y3]

        #4行目の格子点のx,y座標
        Y4 = Y3-3*sqrt(3)/4
        x4=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y4=[Y4   ,   Y4,   Y4,  Y4,  Y4,   Y4,   Y4]

        #5行目の格子点のx,y座標
        Y5 = Y4-3*sqrt(3)/4
        x5=[-3.0,0.0,3.0]
        y5=[Y5, Y5, Y5]

        #6行目の格子点のx,y座標
        Y6 = Y5-3*sqrt(3)/4
        x6=[-21/4,-15/4, -9/4, -3/4, 3/4, 9/4, 15/4]
        y6=[Y6   ,   Y6,   Y6,  Y6,  Y6,   Y6,   Y6]

       #7行目の格子点のx,y座標
        Y7 = Y6-3*sqrt(3)/4
        x7=[-9/2,-3/2, 3/2]
        y7=[ Y7,Y7,Y7]

        X=Float64[]
        Y=Float64[]
        X=copy(x1)
        append!(X,x2)
        append!(X,x3)
        append!(X,x4)
        append!(X,x5)
        append!(X,x6)
        append!(X,x7)
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)

        #2列目のx座標(y座標はYと同じ)
        X_2 = X.+9
        #3列目のx座標(y座標はYと同じ)
        X_3 = X.+18

        #2行目のx,y座標
        X_21 = X .-4.5
        X_22 = X_2 .-4.5
        X_23 = X_3 .-4.5
        Y_2= Y .-(6*3*sqrt(3)/4)

        #3行目のx,y座標
        X_31 = X_21 .-4.5
        X_32 = X_22 .-4.5
        X_33 = X_23 .-4.5
        Y_3 = Y_2 .-(6*3*sqrt(3)/4)

        bond_from=[1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29]
        bond_to=[5,6,7,8,9,10,5,11,6,11,7,12,8,12,9,13,10,13,15,16,17,18,19,20,15,21,16,21,17,22,18,22,19,23,20,23,25,26,27,28,29,30,25,31,26,31,27,32,28,32,29,33,30,33]

        sz=zeros(33)
        open(input_sz_filename) do sz_file
           sz_index=1
           for line in eachline(sz_file)
              s = split(line,",")
              sz[sz_index]=parse(Float64,s[2])    
              sz_index+=1
            end
        end

        szz_rel=zeros(54)
        szz_c_rel = zeros(54)
        #CSVファイルから<S_i^zS_j^z>のデータと<S_i^z>のデータを読み込む
        open(input_szz_filename) do szz_file
            szz_rel_index=1
           for line in eachline(szz_file)
              s = split(line,",")
              szz_rel[szz_rel_index] = parse(Float64, s[3])
                szz_c_rel[szz_rel_index] = parse(Float64, s[4])
                szz_rel_index+=1
            end
        end

        #格子点のplot
         plt_szz=Plots.plot(title=L"$\angle S_i^zS_j^z \rangle$",showaxis=false,margin=0mm)
            title_c = L"$\langle \hat{S_i^z}\hat{S_j^z} \rangle_C$"
         plt_szz_c=Plots.plot(title=title_c,titleposition=:center,showaxis=false,margin=0mm)

        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:red
                m_size=sz[l] * marker_magnification
                plot!(plt_szz,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:blue        
                m_size=abs(sz[l]) * marker_magnification
                plot!(plt_szz,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end

        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:red
                m_size = sz[l] * marker_magnification
                plot!(plt_szz_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:blue        
                m_size = abs(sz[l]) * marker_magnification
                plot!(plt_szz_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end

        #<S_i^zS_j^z>のplot
        for i=1:length(bond_from)
            if szz_rel[i] >= 0
                COLOR=:red
                LW=szz_rel[i]*line_magnification
                plot!(plt_szz, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                COLOR=:blue
                LW=abs(szz_rel[i])*line_magnification
                plot!(plt_szz, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end
        savefig(output_szz_filename)

        #<S_i^zS_j^z> - <S_i^z><S_j^z>のplot
        for i=1:length(bond_from)
            if szz_c_rel[i] >= 0
                COLOR=:red
                LW=abs(szz_c_rel[i])*line_magnification
                plot!(plt_szz_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line,lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                COLOR=:blue
                LW=abs(szz_c_rel[i])*line_magnification
                plot!(plt_szz_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end

        savefig(output_szz_c_filename)
        return plt_szz, plt_szz_c
    end

    function plot_27site_sxx_NNcorellation(input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_sxx_filename)
        #1行目の格子点のx,y座標
        Y1 = 0
        x1=[0.0,3.0,6.0]
        y1=[0.0,0.0,0.0]

        #2行目の格子点のx,y座標
        Y2 = Y1-3*sqrt(3)/4
        x2=[-9/4,-3/4,3/4,9/4,15/4,21/4,27/4]
        y2=[Y2,    Y2,  Y2,Y2,  Y2,  Y2,  Y2]

        #3行目の格子点のx,y座標
        Y3 = Y2-3*sqrt(3)/4
        x3=[-3/2, 3/2,9/2]
        y3=[Y3,Y3,Y3]

        #4行目の格子点のx,y座標
        Y4 = Y3-3*sqrt(3)/4
        x4=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y4=[Y4   ,   Y4,   Y4,  Y4,  Y4,   Y4,   Y4]

        #5行目の格子点のx,y座標
        Y5 = Y4-3*sqrt(3)/4
        x5=[-3.0,0.0,3.0]
        y5=[Y5, Y5, Y5]

        #6行目の格子点のx,y座標
        Y6 = Y5-3*sqrt(3)/4
        x6=[-21/4,-15/4, -9/4, -3/4, 3/4, 9/4, 15/4]
        y6=[Y6   ,   Y6,   Y6,  Y6,  Y6,   Y6,   Y6]

        #7行目の格子点のx,y座標
        Y7 = Y6-3*sqrt(3)/4
        x7=[-9/2,-3/2, 3/2]
        y7=[Y7,Y7,Y7]

        X=Float64[]
        Y=Float64[]
        X=copy(x1)
        append!(X,x2)
        append!(X,x3)
        append!(X,x4)
        append!(X,x5)
        append!(X,x6)
        append!(X,x7)
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)

        #2列目のx座標(y座標はYと同じ)
        X_2 = X.+9
        #3列目のx座標(y座標はYと同じ)
        X_3 = X.+18

        #2行目のx,y座標
        X_21 = X .-4.5
        X_22 = X_2 .-4.5
        X_23 = X_3 .-4.5
        Y_2= Y .-(6*3*sqrt(3)/4)

        #3行目のx,y座標
        X_31 = X_21 .-4.5
        X_32 = X_22 .-4.5
        X_33 = X_23 .-4.5
        Y_3 = Y_2 .-(6*3*sqrt(3)/4)

        bond_from=[1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29]
        bond_to=[5,6,7,8,9,10,5,11,6,11,7,12,8,12,9,13,10,13,15,16,17,18,19,20,15,21,16,21,17,22,18,22,19,23,20,23,25,26,27,28,29,30,25,31,26,31,27,32,28,32,29,33,30,33]

        sz=zeros(33)
        open(input_sz_filename) do sz_file
            sz_index=1
            for line in eachline(sz_file)
                s = split(line,",")
                sz[sz_index]=parse(Float64,s[2])    
                sz_index+=1
            end
        end

        sxx_rel=zeros(54)
        #CSVファイルから<S_i^zS_j^z>のデータと<S_i^z>のデータを読み込む
        open(input_sxx_filename) do sxx_file
            sxx_rel_index=1
            for line in eachline(sxx_file)
                s = split(line,",")
                sxx_rel[sxx_rel_index] = parse(Float64, s[3])
                sxx_rel_index+=1
            end
        end

        #格子点のplot
        title= L"$\langle \hat{S_i^x}\hat{S_j^x} \rangle_C$"
        plt_sxx=Plots.plot(title=title,titleposition=:center,titlevspan=0.2,showaxis=false,margin=0mm)
    
        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:red
                m_size=sz[l] * marker_magnification
                plot!(plt_sxx,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:blue        
                m_size=abs(sz[l]) * marker_magnification
                plot!(plt_sxx,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end

        #<S_ixzS_j^x>のplot
        for i=1:length(bond_from)
            if sxx_rel[i] >= 0
                COLOR=:red
                LW=sxx_rel[i]*line_magnification
                plot!(plt_sxx, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                COLOR=:blue
                LW=abs(sxx_rel[i])*line_magnification
                plot!(plt_sxx, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end
        savefig(output_sxx_filename)

        return plt_sxx
    end

    function plot_36site_szz_NNcorellation(input_szz_filename,input_sz_filename ,marker_magnification,line_magnification,output_szz_filename,output_szz_c_filename)
        Y1=0.0
        x1=[0.0]
        y1=[Y1]
    
        #2行目の格子点のx,y座標
        Y2=Y1-3*sqrt(3)/4
        x2=[-9/4, -3/4, 3/4, 9/4]
        y2=[Y2  ,   Y2,  Y2,  Y2]
    
        #3行目の格子点のx,y座標
        Y3=Y2-3*sqrt(3)/4
        x3=[-9/2, -3/2, 3/2]
        y3=[Y3  ,   Y3,  Y3]
    
        #4行目の格子点のx,y座標
        Y4=Y3-3*sqrt(3)/4
        x4=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3, -27/4+6, -27/4 + (3/2)*5, -27/4 + 9]
        y4=[Y4  ,           Y4,      Y4,            Y4,      Y4,              Y4, Y4]
    
        #5行目の格子点のx,y座標
        Y5=Y4-3*sqrt(3)/4
        x5=[-6,-3,0]
        y5=[Y5,Y5,Y5]
    
        #6行目の格子点のx,y座標
        Y6=Y5-3*sqrt(3)/4
        x6=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3, -27/4+6, -27/4 + (3/2)*5, -27/4 + 9]
        y6=[Y6  ,           Y6,      Y6,            Y6,      Y6,              Y6,        Y6]
    
        #7行目の格子点のx,y座標
        Y7=Y6-3*sqrt(3)/4
        x7=[-9/2, -3/2, 3/2]
        y7=[Y7,Y7,Y7]
    
        #8行目の格子点のx,y座標
        Y8=Y7-3*sqrt(3)/4
        x8=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3, -27/4+6, -27/4 + (3/2)*5, -27/4 + 9]
        y8=[Y8 ,           Y8,      Y8,            Y8,      Y8,              Y8, Y8]
    
        #9行目の格子点のx,y座標
        Y9=Y8-3*sqrt(3)/4
        x9=[-6,-3,0]
        y9=[Y9,Y9,Y9]
    
        #10行目の格子点のx,y座標
        Y10=Y9-3*sqrt(3)/4
        x10=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3]
        y10=[Y10 ,           Y10,      Y10,            Y10]
    
        #11行目の格子点のx,y座標
        Y11=Y10-3*sqrt(3)/4
        x11=[-9/2]
        y11=[Y11]
    
        X=Float64[]
        Y=Float64[]
        X=copy(x1)
        append!(X,x2)
        append!(X,x3)
        append!(X,x4)
        append!(X,x5)
        append!(X,x6)
        append!(X,x7)
        append!(X,x8)
        append!(X,x9)
        append!(X,x10)
        append!(X,x11)
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)
        append!(Y,y8)
        append!(Y,y9)
        append!(Y,y10)
        append!(Y,y11)
    


        #最近接bondの情報を用意する
        bond_from=[1,1,2,2,3,3,4,4,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,16,16,17,17,18,18,19,20,20,21,21,22,22,23,23,24,24,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,36,36,37,37,39,40,40,41,41]
        bond_to=[3,4,3,7,4,7,5,8,8,10,11,12,13,14,15,10,16,11,16,12,17,13,17,14,18,15,18,19,20,21,22,23,24,20,21,26,22,26,23,27,24,27,25,28,28,30,31,32,33,34,35,30,36,31,36,32,37,33,37,34,38,35,38,39,40,41,42,40,41,43,42,43]

        for i=1:length(bond_to)
            println(bond_from[i],",",bond_to[i])
        end

        sz=zeros(43)
        open(input_sz_filename) do sz_file
            sz_index=1
            for line in eachline(sz_file)
                s = split(line,",")
                sz[sz_index]=parse(Float64,s[2])    
                sz_index+=1
            end
        end
  

        szz_rel=zeros(72) #72 = length(bond_from) = length(bond_to)
        szz_c_rel = zeros(72)
        # CSVファイルから<S_i^zS_j^z>のデータと<S_i^z>のデータを読み込む
        open(input_szz_filename) do szz_file
            szz_rel_index=1
            for line in eachline(szz_file)
                s = split(line,",")
                szz_rel[szz_rel_index] = parse(Float64, s[3])
                szz_c_rel[szz_rel_index] = parse(Float64, s[4])
                szz_rel_index+=1
            end
        end
   
        plt=Plots.plot()
        plt_c=Plots.plot()
    
    
        #格子点のplot
        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:red
                m_size=sz[l] * marker_magnification
                plot!(plt,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:blue        
                m_size=abs(sz[l]) * marker_magnification
                plot!(plt,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end
    
        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:red
                m_size = sz[l] * marker_magnification
                plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:blue        
                m_size = abs(sz[l]) * marker_magnification
                plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end
    
        #<S_i^zS_j^z>のplot
        for i=1:length(bond_from)
            if szz_c_rel[i] >= 0
                COLOR=:red
                LW=szz_rel[i]*line_magnification
                plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                COLOR=:blue
                LW=abs(szz_rel[i])*line_magnification
                plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end
        savefig(output_szz_filename)
    
        #<S_i^zS_j^z> - <S_i^z><S_j^z>のplot
        for i=1:length(bond_from)
            if szz_c_rel[i] >= 0
                COLOR=:red
                LW=abs(szz_rel[i])*line_magnification
                plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line,lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                COLOR=:blue
                LW=abs(szz_c_rel[i])*line_magnification
                plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end
    
        savefig(output_szz_c_filename)
    end

    function plot_36site_sxx_NNcorellation(input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_sxx_filename)
        Y1=0.0
        x1=[0.0]
        y1=[Y1]
    
        #2行目の格子点のx,y座標
        Y2=Y1-3*sqrt(3)/4
        x2=[-9/4, -3/4, 3/4, 9/4]
        y2=[Y2  ,   Y2,  Y2,  Y2]
    
        #3行目の格子点のx,y座標
        Y3=Y2-3*sqrt(3)/4
        x3=[-9/2, -3/2, 3/2]
        y3=[Y3  ,   Y3,  Y3]
    
        #4行目の格子点のx,y座標
        Y4=Y3-3*sqrt(3)/4
        x4=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3, -27/4+6, -27/4 + (3/2)*5, -27/4 + 9]
        y4=[Y4  ,           Y4,      Y4,            Y4,      Y4,              Y4, Y4]
    
        #5行目の格子点のx,y座標
        Y5=Y4-3*sqrt(3)/4
        x5=[-6,-3,0]
        y5=[Y5,Y5,Y5]
    
        #6行目の格子点のx,y座標
        Y6=Y5-3*sqrt(3)/4
        x6=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3, -27/4+6, -27/4 + (3/2)*5, -27/4 + 9]
        y6=[Y6  ,           Y6,      Y6,            Y6,      Y6,              Y6,        Y6]
    
        #7行目の格子点のx,y座標
        Y7=Y6-3*sqrt(3)/4
        x7=[-9/2, -3/2, 3/2]
        y7=[Y7,Y7,Y7]
    
        #8行目の格子点のx,y座標
        Y8=Y7-3*sqrt(3)/4
        x8=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3, -27/4+6, -27/4 + (3/2)*5, -27/4 + 9]
        y8=[Y8 ,           Y8,      Y8,            Y8,      Y8,              Y8, Y8]
    
        #9行目の格子点のx,y座標
        Y9=Y8-3*sqrt(3)/4
        x9=[-6,-3,0]
        y9=[Y9,Y9,Y9]
    
        #10行目の格子点のx,y座標
        Y10=Y9-3*sqrt(3)/4
        x10=[-27/4, -27/4+(3/2), -27/4+3, -27/4+(3/2)*3]
        y10=[Y10 ,           Y10,      Y10,            Y10]
    
        #11行目の格子点のx,y座標
        Y11=Y10-3*sqrt(3)/4
        x11=[-9/2]
        y11=[Y11]
    
        X=Float64[]
        Y=Float64[]
        X=copy(x1)
        append!(X,x2)
        append!(X,x3)
        append!(X,x4)
        append!(X,x5)
        append!(X,x6)
        append!(X,x7)
        append!(X,x8)
        append!(X,x9)
        append!(X,x10)
        append!(X,x11)
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)
        append!(Y,y8)
        append!(Y,y9)
        append!(Y,y10)
        append!(Y,y11)
    
        println(length(X))
    
        #最近接bondの情報を用意する
        bond_from=[1,1,2,2,3,3,4,4,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,16,16,17,17,18,18,19,20,20,21,21,22,22,23,23,24,24,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,36,36,37,37,39,40,40,41,41]
        bond_to=[3,4,3,7,4,7,5,8,8,10,11,12,13,14,15,10,16,11,16,12,17,13,17,14,18,15,18,19,20,21,22,23,24,20,21,26,22,26,23,27,24,27,25,28,28,30,31,32,33,34,35,30,36,31,36,32,37,33,37,34,38,35,38,39,40,41,42,40,41,43,42,43]
   
        sz=zeros(43)
        open(input_sz_filename) do sz_file
            sz_index=1
            for line in eachline(sz_file)
                s = split(line,",")
                sz[sz_index]=parse(Float64,s[2])    
                sz_index+=1
            end
        end
  

        sxx_rel=zeros(72) #72 = length(bond_from) = length(bond_to)
        # CSVファイルから<S_i^zS_j^z>のデータと<S_i^z>のデータを読み込む
        open(input_sxx_filename) do sxx_file
            sxx_rel_index=1
            for line in eachline(sxx_file)
                s = split(line,",")
                sxx_rel[sxx_rel_index] = parse(Float64, s[3])
                sxx_rel_index+=1
            end
        end

        plt=Plots.plot() 
    
        #格子点のplot
        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:red
                m_size=sz[l] * marker_magnification
            plot!(plt,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:blue        
                m_size=abs(sz[l]) * marker_magnification
                plot!(plt,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end
    
        #<S_i^zS_j^z>のplot
        for i=1:length(bond_from)
        if sxx_rel[i] >= 0
            COLOR=:red
            LW=sxx_rel[i]*line_magnification
            plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
        else
            COLOR=:blue
            LW=abs(sxx_rel[i])*line_magnification
            plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
        end
    end
    savefig(output_sxx_filename)
end

export plot_27site_MH_curve_with_Highrigt
export plot_27site_szz_NNcorellation
export plot_27site_sxx_NNcorellation
export plot_36site_szz_NNcorellation
export plot_36site_sxx_NNcorellation
end
