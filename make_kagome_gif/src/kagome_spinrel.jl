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

        h[15] = 10
        M[15] = M[14]

        plt = plot(h,M,line=:steppost,xlabel=L"h", ylabel=L"M/M_{sat}",xlims=(0,3.5),ylims=(0,1),title="Magnetization curve",titleposition=:center,legend=false, titlegap=2)
        if num_of_M < 27
            plot!(plt, [h[num_of_M-13],h[num_of_M-12]],[M[num_of_M-13],M[num_of_M-13]],line=:steppost, lw=3,linecolor=:red)
        elseif num_of_M == 27
            plot!(plt, [h[num_of_M-13],3.5],[M[num_of_M-13],M[num_of_M-13]],line=:steppost,lw=3, linecolor=:red, titlegap=2)
        end
        return plt
    end

    function plot_36site_MH_curve_with_Highrigt(MHdata_file, num_of_M)
        h=zeros(20)
        M=zeros(20)
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

        h[20] = 10
        M[20] = M[19]

        plt = plot(h,M,line=:steppost,xlabel=L"h", ylabel=L"M/M_{sat}",xlims=(0,3.5),ylims=(0,1),title="Magnetization curve",titleposition=:center,legend=false, titlegap=2)
        if num_of_M < 36
            plot!(plt, [h[num_of_M-17],h[num_of_M-16]],[M[num_of_M-17],M[num_of_M-17]],line=:steppost, lw=3,linecolor=:red, titlegap=2)
        elseif num_of_M == 36
            plot!(plt, [h[num_of_M-17],3.5],[M[num_of_M-17],M[num_of_M-17]],line=:steppost,lw=3, linecolor=:red, titlegap=2)
        end
        return plt
    end

    function plot_27site_kagome_lattice(J_red, J_green, J_blue)
        #1行目の格子点のx,y座標
        Y1 = 0.0
        # x1 = [0.0, 3/2, 3, 9/2]
        # y1 = [0.0,0.0,0.0,0.0]
        # Y1 = 0
        x1=[0.0,3.0,6.0]
        y1=[0.0,0.0,0.0]

        #2行目の格子点のx,y座標
        Y2 = Y1-3*sqrt(3)/4
        # x2 = [-3/4,9/4, 21/4]
        # y2 = [Y2,Y2,Y2]
        x2=[-9/4,-3/4,3/4,9/4,15/4,21/4,27/4]
        y2=[Y2,    Y2,  Y2,Y2,  Y2,  Y2,  Y2]

        #3行目の格子点のx,y座標
        Y3 = Y2-3*sqrt(3)/4
        # x3 = [-3/2, 0, 3/2, 3, 9/2, 6]
        # y3 = [Y3,Y3,Y3,Y3,Y3,Y3]
        x3=[-3/2, 3/2,9/2]
        y3=[Y3,Y3,Y3]

        #4行目の格子点のx,y座標
        Y4 = Y3-3*sqrt(3)/4
        # x4 = [-9/4, 3/4, 15/4, 27/4]
        # y4 = [Y4   ,   Y4,   Y4,  Y4]
        x4=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y4=[Y4   ,   Y4,   Y4,  Y4,  Y4,   Y4,   Y4]

        #5行目の格子点のx,y座標
        Y5 = Y4-3*sqrt(3)/4
        # x5 = [-3/2, 0, 3/2, 3, 9/2, 6]
        # y5 = [Y5  ,Y5,  Y5, Y5, Y5, Y5]
        x5=[-3.0,0.0,3.0]
        y5=[Y5, Y5, Y5]

        #6行目の格子点のx,y座標
        Y6 = Y5-3*sqrt(3)/4
        # x6 = [-3/4,9/4, 21/4]
        # y6 = [Y6, Y6, Y6]
        x6=[-21/4,-15/4, -9/4, -3/4, 3/4, 9/4, 15/4]
        y6=[Y6   ,   Y6,   Y6,  Y6,  Y6,   Y6,   Y6]

       #7行目の格子点のx,y座標
        Y7 = Y6-3*sqrt(3)/4
        # x7 = [0.0, 3/2, 3, 9/2]
        # y7 = [Y7,Y7,Y7,Y7]
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
        
        #4列目のx座標
        X_4 = X.+27
        
        #2行目のx,y座標
        X_21 = X .-4.5
        X_22 = X_2 .-4.5
        X_23 = X_3 .-4.5
        X_24 = X_4 .-4.5
        Y_2= Y .-(6*3*sqrt(3)/4)

        #3行目のx,y座標
        X_31 = X_21 .-4.5
        X_32 = X_22 .-4.5
        X_33 = X_23 .-4.5
        Y_3 = Y_2 .-(6*3*sqrt(3)/4)

        bond_from=[1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29]
        bond_to=[5,6,7,8,9,10,5,11,6,11,7,12,8,12,9,13,10,13,15,16,17,18,19,20,15,21,16,21,17,22,18,22,19,23,20,23,25,26,27,28,29,30,25,31,26,31,27,32,28,32,29,33,30,33]
        bond_strength=[J_red, J_blue, J_blue, J_green, J_green, J_red,J_red,J_green,J_green,J_blue,J_green,J_blue,J_red,J_red,J_blue,J_red,
                    J_blue,J_green,J_blue,J_green,J_green,J_red,J_red,J_blue,J_green,J_blue,J_red,J_red,J_blue,J_red,J_blue,J_green,J_red,J_green,J_green,
                    J_blue,J_green,J_red,J_red,J_blue,J_blue,J_green,J_blue,J_red,J_blue,J_green,J_red,J_green,J_green,J_blue,J_green,J_blue,J_red,J_red]

        
        sz=ones(33)

        plt=Plots.plot(title=L"$\mathrm{Kagome\ Lattice,\ 27site}$",showaxis=false)
        
        #<S_i^zS_j^z>のplot
        for i=1:length(bond_from)
            if bond_strength[i] == J_red
                COLOR=:red
             elseif bond_strength[i] == J_green
                COLOR=:green
            elseif bond_strength[i] == J_blue
                COLOR=:blue
            end

            LW=bond_strength[i]*5
            plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
        end
    
        #格子点のplot
        for l=1:length(sz)
            m_color=:black
            m_size=sz[l] * 5
            plot!(plt,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)           
            plot!(plt,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        end

        return plt
        
    end


    function plot_36site_kagome_lattice(J_red,J_green,J_blue)
        Y1=0.0
        x1=[0.0, 3.0]
        y1=[Y1,Y1]
    
        #2行目の格子点のx,y座標
        Y2=Y1-3*sqrt(3)/4
        x2=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y2=[Y2  ,   Y2,  Y2,  Y2,   Y2,   Y2]
    
        #3行目の格子点のx,y座標
        Y3=Y2-3*sqrt(3)/4
        x3=[-3/2, 3/2, 9/2]
        y3=[Y3  ,   Y3, Y3]
    
        #4行目の格子点のx,y座標
        Y4=Y3-3*sqrt(3)/4
        x4=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y4=[Y4  ,    Y4,   Y4,  Y4,  Y4,   Y4,  Y4,    Y4]
    
        #5行目の格子点のx,y座標
        Y5=Y4-3*sqrt(3)/4
        x5=[-3,0, 3, 6]
        y5=[Y5,Y5,Y5,Y5]
    
        #6行目の格子点のx,y座標
        Y6=Y5-3*sqrt(3)/4
        x6=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y6=[Y6  ,    Y6,   Y6,  Y6,  Y6,   Y6,   Y6,   Y6]
    
        #7行目の格子点のx,y座標
        Y7=Y6-3*sqrt(3)/4
        x7=[-3/2, 3/2, 9/2]
        y7=[Y7,Y7,Y7]
    
        #8行目の格子点のx,y座標
        Y8=Y7-3*sqrt(3)/4
        x8=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y8=[Y8 ,    Y8,  Y8,  Y8,   Y8,   Y8]
    
        #9行目の格子点のx,y座標
        Y9=Y8-3*sqrt(3)/4
        x9=[0.0, 3.0]
        y9=[Y9,Y9]

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
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)
        append!(Y,y8)
        append!(Y,y9)

        #2列目のx座標(y座標はYと同じ)
        X_2=X.+9
        Y_2=Y.+3*sqrt(3)
        #3列目のx座標
        X_3 = X.+9
        Y_3 = Y.-3*sqrt(3)
        X_4 = X
        Y_4 = Y.-6*sqrt(3)
        X_5 = X_3.+9
        Y_5 = Y_3.+3*sqrt(3)
        X_6 = X_3.+9
        Y_6 = Y_3.-3*sqrt(3)
        X_7 = X_4.+9
        Y_7 = Y_4.-3*sqrt(3)

        bond_from=[1,1 ,2,2,3,3,4,4 ,5  ,5  ,6 ,6 ,7 ,7,8 ,9  ,9   ,10 ,10 ,11,11  ,12    ,12
        ,13        ,13        ,14        ,14        ,15        ,15        ,16        ,16
        ,17        ,17        ,18        ,18        ,19        ,20        ,20        ,21
        ,21        ,22        ,22        ,23        ,23        ,24        ,25        ,25
        ,26        ,26        ,27        ,27        ,28        ,28        ,29        ,29
        ,30        ,30        ,32        ,32        ,33        ,33        ,34        ,34
        ,35        ,36        ,36        ,37        ,37        ,38        ,38        ,39        ,39        ]
        bond_to=[4        ,5        ,6        ,7        ,4        ,9        ,5        ,9
        ,6        ,10        ,7        ,10        ,8        ,11        ,11        ,13
        ,14        ,15        ,16        ,17        ,18        ,13        ,20        ,14
        ,20        ,15        ,21        ,16        ,21       , 17        ,22        ,18
        ,22        ,19        ,23        ,23        ,24        ,25        ,26        ,27
        ,28        ,29        ,30        ,31        ,25        ,26        ,32        ,27
        ,32        ,28        ,33        ,29        ,33        ,30        ,34        ,31
        ,34        ,35        ,36        ,37        ,38        ,39        ,40        ,36
        ,37        ,41        ,38        ,41        ,39        ,42       ,40        ,42        ]

        bond_strength=[J_green        ,J_red        ,J_red        ,J_blue        ,J_blue        ,J_red        ,J_blue        ,J_green
        ,J_red        ,J_green        ,J_green        ,J_blue        ,J_green        ,J_blue        ,J_red        ,J_red
        ,J_blue        ,J_blue        ,J_green        ,J_green        ,J_red        ,J_red        ,J_green        ,J_green
        ,J_blue        ,J_green        ,J_blue        ,J_red        ,J_red        ,J_blue        ,J_red        ,J_blue
        ,J_green        ,J_red        ,J_green        ,J_blue        ,J_blue        ,J_green        ,J_green        ,J_red
        ,J_red        ,J_blue        ,J_blue        ,J_green        ,J_red        ,J_blue        ,J_red        ,J_blue
        ,J_green        ,J_red        ,J_green        ,J_green        ,J_blue        ,J_green        ,J_blue        ,J_red
        ,J_red        ,J_red        ,J_blue        ,J_blue        ,J_green        ,J_green        ,J_red        ,J_green
        ,J_green        ,J_blue        ,J_red        ,J_red        ,J_blue        ,J_red        ,J_blue        ,J_green        ]
    
        sz = ones(Int,42)
        COLOR=:black
        plt=Plots.plot(title=L"$\mathrm{Kagome Lattice, 36site}$",showaxis=false)
        #bondのplot
        for i=1:length(bond_from)
            if bond_strength[i] == J_red
                COLOR=:red
             elseif bond_strength[i] == J_green
                COLOR=:green
            elseif bond_strength[i] == J_blue
                COLOR=:blue
            end

            LW=bond_strength[i]*5
            plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
        end

    
        #格子点のplot
         for l=1:length(sz)
            m_color=:black
            m_size=sz[l] * 5
            plot!(plt,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal) 
            plot!(plt,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)     
            plot!(plt,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)        
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

        # # #格子点のplot
         plt_szz=Plots.plot(title=L"$\langle S_i^zS_j^z \rangle$",showaxis=false)
            title_c = L"$\langle \hat{S_i^z}\hat{S_j^z} \rangle_C$"
         plt_szz_c=Plots.plot(title=title_c,titleposition=:center,showaxis=false)

        #          #<S_i^zS_j^z>のplot
        for i=1:length(bond_from)
            if szz_rel[i] >= 0
                COLOR=:violetred2
                LW=szz_rel[i]*line_magnification
                plot!(plt_szz, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                COLOR=:dodgerblue1
                LW=abs(szz_rel[i])*line_magnification
                plot!(plt_szz, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end


        savefig(output_szz_c_filename)

        # # #<S_i^zS_j^z> - <S_i^z><S_j^z>のplot
        for i=1:length(bond_from)
            if szz_c_rel[i] >= 0
                COLOR=:violetred2
                LW=abs(szz_c_rel[i])*line_magnification
                plot!(plt_szz_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line,lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                COLOR=:dodgerblue1
                LW=abs(szz_c_rel[i])*line_magnification
                plot!(plt_szz_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end

        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:violetred2
                m_size=sz[l] * marker_magnification
                plot!(plt_szz,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:dodgerblue1        
                m_size=abs(sz[l]) * marker_magnification
                plot!(plt_szz,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end

        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:violetred2
                m_size = sz[l] * marker_magnification
                plot!(plt_szz_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:dodgerblue1        
                m_size = abs(sz[l]) * marker_magnification
                plot!(plt_szz_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_szz_c,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_szz_c,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
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
        title= L"$\langle \hat{S_i^x}\hat{S_j^x} \rangle$"
        plt_sxx=Plots.plot(title=title,titleposition=:center,titlevspan=0.2,showaxis=false)

                #<S_ixzS_j^x>のplot
                for i=1:length(bond_from)
                    if sxx_rel[i] >= 0
                        COLOR=:violetred2
                        LW=sxx_rel[i]*line_magnification
                        plot!(plt_sxx, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        plot!(plt_sxx, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        plot!(plt_sxx, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        plot!(plt_sxx, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                    else
                        COLOR=:dodgerblue1
                        LW=abs(sxx_rel[i])*line_magnification
                        plot!(plt_sxx, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        plot!(plt_sxx, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        plot!(plt_sxx, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        plot!(plt_sxx, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_23[bond_from[i]],X_23[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_31[bond_from[i]],X_31[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_32[bond_from[i]],X_32[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                        # plot!(plt_sxx, [X_33[bond_from[i]],X_33[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                        # framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                    end
                end

        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:violetred2
                m_size=sz[l] * marker_magnification
                plot!(plt_sxx,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            else
                m_color=:dodgerblue1        
                m_size=abs(sz[l]) * marker_magnification
                plot!(plt_sxx,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_3[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_sxx,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_23[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_31[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_32[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                # plot!(plt_sxx,[X_33[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
            end
        end


        savefig(output_sxx_filename)

        return plt_sxx
    end

    function plot_36site_szz_NNcorellation(input_szz_filename,input_sz_filename ,marker_magnification,line_magnification,output_szz_filename,output_szz_c_filename)
        Y1=0.0
        x1=[0.0, 3.0]
        y1=[Y1,Y1]
    
        #2行目の格子点のx,y座標
        Y2=Y1-3*sqrt(3)/4
        x2=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y2=[Y2  ,   Y2,  Y2,  Y2,   Y2,   Y2]
    
        #3行目の格子点のx,y座標
        Y3=Y2-3*sqrt(3)/4
        x3=[-3/2, 3/2, 9/2]
        y3=[Y3  ,   Y3, Y3]
    
        #4行目の格子点のx,y座標
        Y4=Y3-3*sqrt(3)/4
        x4=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y4=[Y4  ,    Y4,   Y4,  Y4,  Y4,   Y4,  Y4,    Y4]
    
        #5行目の格子点のx,y座標
        Y5=Y4-3*sqrt(3)/4
        x5=[-3,0, 3, 6]
        y5=[Y5,Y5,Y5,Y5]
    
        #6行目の格子点のx,y座標
        Y6=Y5-3*sqrt(3)/4
        x6=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y6=[Y6  ,    Y6,   Y6,  Y6,  Y6,   Y6,   Y6,   Y6]
    
        #7行目の格子点のx,y座標
        Y7=Y6-3*sqrt(3)/4
        x7=[-3/2, 3/2, 9/2]
        y7=[Y7,Y7,Y7]
    
        #8行目の格子点のx,y座標
        Y8=Y7-3*sqrt(3)/4
        x8=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y8=[Y8 ,    Y8,  Y8,  Y8,   Y8,   Y8]
    
        #9行目の格子点のx,y座標
        Y9=Y8-3*sqrt(3)/4
        x9=[0.0, 3.0]
        y9=[Y9,Y9]

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
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)
        append!(Y,y8)
        append!(Y,y9)

        #2列目のx座標(y座標はYと同じ)
        X_2=X.+9
        Y_2=Y.+3*sqrt(3)
        #3列目のx座標
        X_3 = X.+9
        Y_3 = Y.-3*sqrt(3)
        X_4 = X
        Y_4 = Y.-6*sqrt(3)
        X_5 = X_3.+9
        Y_5 = Y_3.+3*sqrt(3)
        X_6 = X_3.+9
        Y_6 = Y_3.-3*sqrt(3)
        X_7 = X_4.+9
        Y_7 = Y_4.-3*sqrt(3)

        bond_from=[1,1 ,2,2,3,3,4,4 ,5  ,5  ,6 ,6 ,7 ,7,8 ,9  ,9   ,10 ,10 ,11,11  ,12    ,12
        ,13        ,13        ,14        ,14        ,15        ,15        ,16        ,16
        ,17        ,17        ,18        ,18        ,19        ,20        ,20        ,21
        ,21        ,22        ,22        ,23        ,23        ,24        ,25        ,25
        ,26        ,26        ,27        ,27        ,28        ,28        ,29        ,29
        ,30        ,30        ,32        ,32        ,33        ,33        ,34        ,34
        ,35        ,36        ,36        ,37        ,37        ,38        ,38        ,39        ,39        ]
        bond_to=[4        ,5        ,6        ,7        ,4        ,9        ,5        ,9
        ,6        ,10        ,7        ,10        ,8        ,11        ,11        ,13
        ,14        ,15        ,16        ,17        ,18        ,13        ,20        ,14
        ,20        ,15        ,21        ,16        ,21       , 17        ,22        ,18
        ,22        ,19        ,23        ,23        ,24        ,25        ,26        ,27
        ,28        ,29        ,30        ,31        ,25        ,26        ,32        ,27
        ,32        ,28        ,33        ,29        ,33        ,30        ,34        ,31
        ,34        ,35        ,36        ,37        ,38        ,39        ,40        ,36
        ,37        ,41        ,38        ,41        ,39        ,42       ,40        ,42        ]

        
        sz=zeros(42)
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
   
        plt=Plots.plot(title=L"$\angle S_i^zS_j^z \rangle$",showaxis=false)
        title_c = L"$\langle \hat{S_i^z}\hat{S_j^z} \rangle_C$"
         plt_c=Plots.plot(title=title_c,titleposition=:center,showaxis=false)
    
            #<S_i^zS_j^z>のplot : 必要であれば以下をコメントアウトして使用する
            # for i=1:length(bond_from)
            #     if szz_rel[i] >= 0
            #         COLOR=:violetred2
            #         LW=szz_rel[i]*line_magnification
            #         plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #     else
            #         COLOR=:dodgerblue1
            #         LW=abs(szz_rel[i])*line_magnification
            #         plot!(plt, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #         plot!(plt, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            #         framestyle=:none,aspect_ratio=:equal)
            #     end
            # end
            # savefig(output_szz_filename)
        
            #<S_i^zS_j^z> - <S_i^z><S_j^z>のplot
            for i=1:length(bond_from)
                if szz_c_rel[i] >= 0
                    COLOR=:violetred2
                    LW=abs(szz_c_rel[i])*line_magnification
                    plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line,lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                else
                    COLOR=:dodgerblue1
                    LW=abs(szz_c_rel[i])*line_magnification
                    plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                    plot!(plt_c, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                end
            end

    
        #格子点のplot
        #<S_i^zS_j^z>のplotが必要であればコメントアウトして使用する
        # for l=1:length(sz)
        #     if sz[l] >= 0
        #         m_color=:violetred2
        #         m_size=sz[l] * marker_magnification
        #         plot!(plt,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)              
        #         plot!(plt,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #     else
        #         m_color=:dodgerblue1        
        #         m_size=abs(sz[l]) * marker_magnification
        #         plot!(plt,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)              
        #         plot!(plt,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #         plot!(plt,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        #     end
        # end
    
        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:violetred2
                m_size = sz[l] * marker_magnification
                plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            else
                m_color=:dodgerblue1        
                m_size = abs(sz[l]) * marker_magnification
                plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,background_color=:transparent,aspect_ratio=:equal)
                plot!(plt_c,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_c,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            end
        end
    
        savefig(output_szz_c_filename)
        return plt, plt_c
    end

    function plot_36site_sxx_NNcorellation(input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_sxx_filename)
        Y1=0.0
        x1=[0.0, 3.0]
        y1=[Y1,Y1]
    
        #2行目の格子点のx,y座標
        Y2=Y1-3*sqrt(3)/4
        x2=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y2=[Y2  ,   Y2,  Y2,  Y2,   Y2,   Y2]
    
        #3行目の格子点のx,y座標
        Y3=Y2-3*sqrt(3)/4
        x3=[-3/2, 3/2, 9/2]
        y3=[Y3  ,   Y3, Y3]
    
        #4行目の格子点のx,y座標
        Y4=Y3-3*sqrt(3)/4
        x4=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y4=[Y4  ,    Y4,   Y4,  Y4,  Y4,   Y4,  Y4,    Y4]
    
        #5行目の格子点のx,y座標
        Y5=Y4-3*sqrt(3)/4
        x5=[-3,0, 3, 6]
        y5=[Y5,Y5,Y5,Y5]
    
        #6行目の格子点のx,y座標
        Y6=Y5-3*sqrt(3)/4
        x6=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y6=[Y6  ,    Y6,   Y6,  Y6,  Y6,   Y6,   Y6,   Y6]
    
        #7行目の格子点のx,y座標
        Y7=Y6-3*sqrt(3)/4
        x7=[-3/2, 3/2, 9/2]
        y7=[Y7,Y7,Y7]
    
        #8行目の格子点のx,y座標
        Y8=Y7-3*sqrt(3)/4
        x8=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y8=[Y8 ,    Y8,  Y8,  Y8,   Y8,   Y8]
    
        #9行目の格子点のx,y座標
        Y9=Y8-3*sqrt(3)/4
        x9=[0.0, 3.0]
        y9=[Y9,Y9]

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
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)
        append!(Y,y8)
        append!(Y,y9)

        #2列目のx座標(y座標はYと同じ)
        X_2=X.+9
        Y_2=Y.+3*sqrt(3)
        #3列目のx座標
        X_3 = X.+9
        Y_3 = Y.-3*sqrt(3)
        X_4 = X
        Y_4 = Y.-6*sqrt(3)
        X_5 = X_3.+9
        Y_5 = Y_3.+3*sqrt(3)
        X_6 = X_3.+9
        Y_6 = Y_3.-3*sqrt(3)
        X_7 = X_4.+9
        Y_7 = Y_4.-3*sqrt(3)

        bond_from=[1,1 ,2,2,3,3,4,4 ,5  ,5  ,6 ,6 ,7 ,7,8 ,9  ,9   ,10 ,10 ,11,11  ,12    ,12
        ,13        ,13        ,14        ,14        ,15        ,15        ,16        ,16
        ,17        ,17        ,18        ,18        ,19        ,20        ,20        ,21
        ,21        ,22        ,22        ,23        ,23        ,24        ,25        ,25
        ,26        ,26        ,27        ,27        ,28        ,28        ,29        ,29
        ,30        ,30        ,32        ,32        ,33        ,33        ,34        ,34
        ,35        ,36        ,36        ,37        ,37        ,38        ,38        ,39        ,39        ]
        bond_to=[4        ,5        ,6        ,7        ,4        ,9        ,5        ,9
        ,6        ,10        ,7        ,10        ,8        ,11        ,11        ,13
        ,14        ,15        ,16        ,17        ,18        ,13        ,20        ,14
        ,20        ,15        ,21        ,16        ,21       , 17        ,22        ,18
        ,22        ,19        ,23        ,23        ,24        ,25        ,26        ,27
        ,28        ,29        ,30        ,31        ,25        ,26        ,32        ,27
        ,32        ,28        ,33        ,29        ,33        ,30        ,34        ,31
        ,34        ,35        ,36        ,37        ,38        ,39        ,40        ,36
        ,37        ,41        ,38        ,41        ,39        ,42       ,40        ,42        ]


        sz=zeros(42)
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

        title= L"$\langle \hat{S_i^x}\hat{S_j^x} \rangle$"
        plt_sxx=Plots.plot(title=title,titleposition=:center,showaxis=false)
    
        #<S_i^xS_j^x>のplot
        for i=1:length(bond_from)
            if sxx_rel[i] >= 0
                COLOR=:violetred2
                LW=sxx_rel[i]*line_magnification
                plot!(plt_sxx, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
            else
                COLOR=:dodgerblue1
                LW=abs(sxx_rel[i])*line_magnification
                plot!(plt_sxx, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                    framestyle=:none,aspect_ratio=:equal)
            end
        end

        #格子点のplot
        for l=1:length(sz)
            if sz[l] >= 0
                m_color=:violetred2
                m_size=sz[l] * marker_magnification
                plot!(plt_sxx,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)              
                plot!(plt_sxx,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            else
                m_color=:dodgerblue1        
                m_size=abs(sz[l]) * marker_magnification
                plot!(plt_sxx,[X[l]],[Y[l]],label="",st=:scatter,  ms=m_size,msw=0, mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)              
                plot!(plt_sxx,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
                plot!(plt_sxx,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            end
        end

    savefig(output_sxx_filename)

    return plt_sxx
end

function plot_27site_NN_spin_rel(input_szz_filename,input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_spin_rel_filename)
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

    #<S_i^z S_j^z>_cに関するデータのinput
    szz_c_rel = zeros(54)
    #CSVファイルから<S_i^zS_j^z>のデータと<S_i^z>のデータを読み込む
    open(input_szz_filename) do szz_file
        szz_rel_index=1
       for line in eachline(szz_file)
          s = split(line,",")
          szz_c_rel[szz_rel_index] = parse(Float64, s[4])
          szz_rel_index+=1
        end
    end

    #<S_i^x S_j^x>_xに関するデータのinput
    sxx_rel=zeros(54) #54 = length(bond_from) = length(bond_to)
    # CSVファイルから<S_i^zS_j^z>のデータと<S_i^z>のデータを読み込む
    open(input_sxx_filename) do sxx_file
        sxx_rel_index=1
        for line in eachline(sxx_file)
            s = split(line,",")
            sxx_rel[sxx_rel_index] = parse(Float64, s[3])
            sxx_rel_index+=1
        end
    end

    #<S_i S_j>_C = <S_i^z S_j^z>_c - <S_i^x S_j^x>_xの計算
    for i=1:length(szz_c_rel)
        szz_c_rel[i] = szz_c_rel[i] + 2*sxx_rel[i]
    end

    #格子点のplot
    title_c = L"spin_rel"
     plt_c=Plots.plot(title=title_c,titleposition=:center,showaxis=false)

    #<S_iS_j>_Cのplot    
    for i=1:length(bond_from)
        if szz_c_rel[i] >= 0
            COLOR=:violetred2
            LW=abs(szz_c_rel[i])*line_magnification
            plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line,lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
        else
            COLOR=:dodgerblue1
            LW=abs(szz_c_rel[i])*line_magnification
            plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_21[bond_from[i]],X_21[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_22[bond_from[i]],X_22[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
            framestyle=:none,aspect_ratio=:equal)
        end
    end

    for l=1:length(sz)
        if sz[l] >= 0
            m_color=:violetred2
            m_size = sz[l] * marker_magnification
            plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        else
            m_color=:dodgerblue1        
            m_size = abs(sz[l]) * marker_magnification
            plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_2[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_21[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_22[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        end
    end
    savefig(output_spin_rel_filename)
    return plt_c
end

#<S_i S_j>_C = <S_i^z S_j^z> - <S_i^z> - <S_j^z> + 2*<S_i^x S_j^x>
function plot_36site_NN_spin_rel(input_szz_filename,input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_spin_rel_filename)
    Y1=0.0
        x1=[0.0, 3.0]
        y1=[Y1,Y1]
    
        #2行目の格子点のx,y座標
        Y2=Y1-3*sqrt(3)/4
        x2=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y2=[Y2  ,   Y2,  Y2,  Y2,   Y2,   Y2]
    
        #3行目の格子点のx,y座標
        Y3=Y2-3*sqrt(3)/4
        x3=[-3/2, 3/2, 9/2]
        y3=[Y3  ,   Y3, Y3]
    
        #4行目の格子点のx,y座標
        Y4=Y3-3*sqrt(3)/4
        x4=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y4=[Y4  ,    Y4,   Y4,  Y4,  Y4,   Y4,  Y4,    Y4]
    
        #5行目の格子点のx,y座標
        Y5=Y4-3*sqrt(3)/4
        x5=[-3,0, 3, 6]
        y5=[Y5,Y5,Y5,Y5]
    
        #6行目の格子点のx,y座標
        Y6=Y5-3*sqrt(3)/4
        x6=[-15/4, -9/4, -3/4, 3/4, 9/4, 15/4, 21/4, 27/4]
        y6=[Y6  ,    Y6,   Y6,  Y6,  Y6,   Y6,   Y6,   Y6]
    
        #7行目の格子点のx,y座標
        Y7=Y6-3*sqrt(3)/4
        x7=[-3/2, 3/2, 9/2]
        y7=[Y7,Y7,Y7]
    
        #8行目の格子点のx,y座標
        Y8=Y7-3*sqrt(3)/4
        x8=[-9/4, -3/4, 3/4, 9/4, 15/4, 21/4]
        y8=[Y8 ,    Y8,  Y8,  Y8,   Y8,   Y8]
    
        #9行目の格子点のx,y座標
        Y9=Y8-3*sqrt(3)/4
        x9=[0.0, 3.0]
        y9=[Y9,Y9]

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
        Y=copy(y1)
        append!(Y,y2)
        append!(Y,y3)
        append!(Y,y4)
        append!(Y,y5)
        append!(Y,y6)
        append!(Y,y7)
        append!(Y,y8)
        append!(Y,y9)

        #2列目のx座標(y座標はYと同じ)
        X_2=X.+9
        Y_2=Y.+3*sqrt(3)
        #3列目のx座標
        X_3 = X.+9
        Y_3 = Y.-3*sqrt(3)
        X_4 = X
        Y_4 = Y.-6*sqrt(3)
        X_5 = X_3.+9
        Y_5 = Y_3.+3*sqrt(3)
        X_6 = X_3.+9
        Y_6 = Y_3.-3*sqrt(3)
        X_7 = X_4.+9
        Y_7 = Y_4.-3*sqrt(3)

        bond_from=[1,1 ,2,2,3,3,4,4 ,5  ,5  ,6 ,6 ,7 ,7,8 ,9  ,9   ,10 ,10 ,11,11  ,12    ,12
        ,13        ,13        ,14        ,14        ,15        ,15        ,16        ,16
        ,17        ,17        ,18        ,18        ,19        ,20        ,20        ,21
        ,21        ,22        ,22        ,23        ,23        ,24        ,25        ,25
        ,26        ,26        ,27        ,27        ,28        ,28        ,29        ,29
        ,30        ,30        ,32        ,32        ,33        ,33        ,34        ,34
        ,35        ,36        ,36        ,37        ,37        ,38        ,38        ,39        ,39        ]
        bond_to=[4        ,5        ,6        ,7        ,4        ,9        ,5        ,9
        ,6        ,10        ,7        ,10        ,8        ,11        ,11        ,13
        ,14        ,15        ,16        ,17        ,18        ,13        ,20        ,14
        ,20        ,15        ,21        ,16        ,21       , 17        ,22        ,18
        ,22        ,19        ,23        ,23        ,24        ,25        ,26        ,27
        ,28        ,29        ,30        ,31        ,25        ,26        ,32        ,27
        ,32        ,28        ,33        ,29        ,33        ,30        ,34        ,31
        ,34        ,35        ,36        ,37        ,38        ,39        ,40        ,36
        ,37        ,41        ,38        ,41        ,39        ,42       ,40        ,42        ]

    sz=zeros(42)
    open(input_sz_filename) do sz_file
        sz_index=1
        for line in eachline(sz_file)
            s = split(line,",")
            sz[sz_index]=parse(Float64,s[2])    
            sz_index+=1
        end
    end

    #<S_i^z S_j^z>_cに関するデータのinput
    szz_c_rel = zeros(72)  #72 = length(bond_from) = length(bond_to)
    # CSVファイルから<S_i^zS_j^z>のデータと<S_i^z>のデータを読み込む
    open(input_szz_filename) do szz_file
        szz_rel_index=1
        for line in eachline(szz_file)
            s = split(line,",")
            szz_c_rel[szz_rel_index] = parse(Float64, s[4])
            szz_rel_index+=1
        end
    end

    #<S_i^x S_j^x>_xに関するデータのinput
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

    #<S_i S_j>_C = <S_i^z S_j^z>_c + <S_i^x S_j^x>+ <S_i^y S_j^y>の計算
    for i=1:length(szz_c_rel)
        szz_c_rel[i] = szz_c_rel[i] + 2*sxx_rel[i]
    end


    title_c = L"$\langle \hat{\bm{S}_i}\hat{\bm{S}_j} \rangle_C$"
    plt_c=Plots.plot(title=title_c,titleposition=:center,showaxis=false)

    #<S_i S_j>_Cのplot
    for i=1:length(bond_from)
        if szz_c_rel[i] >= 0
            COLOR=:violetred2
            LW=abs(szz_c_rel[i])*line_magnification
            plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line,lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
        else
            COLOR=:dodgerblue1
            LW=abs(szz_c_rel[i])*line_magnification
            plot!(plt_c, [X[bond_from[i]],X[bond_to[i]]],[Y[bond_from[i]],Y[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_2[bond_from[i]],X_2[bond_to[i]]],[Y_2[bond_from[i]],Y_2[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_3[bond_from[i]],X_3[bond_to[i]]],[Y_3[bond_from[i]],Y_3[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_4[bond_from[i]],X_4[bond_to[i]]],[Y_4[bond_from[i]],Y_4[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_5[bond_from[i]],X_5[bond_to[i]]],[Y_5[bond_from[i]],Y_5[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_6[bond_from[i]],X_6[bond_to[i]]],[Y_6[bond_from[i]],Y_6[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c, [X_7[bond_from[i]],X_7[bond_to[i]]],[Y_7[bond_from[i]],Y_7[bond_to[i]]],st=:line, lw=LW,lc=COLOR,label="",
                framestyle=:none,aspect_ratio=:equal)
        end
    end

    #格子点のplot
    for l=1:length(sz)
        if sz[l] >= 0
            m_color=:violetred2
            m_size = sz[l] * marker_magnification
            plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        else
            m_color=:dodgerblue1        
            m_size = abs(sz[l]) * marker_magnification
            plot!(plt_c,[X[l]],[Y[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_2[l]],[Y_2[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_3[l]],[Y_3[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_4[l]],[Y_4[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_5[l]],[Y_5[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_6[l]],[Y_6[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
            plot!(plt_c,[X_7[l]],[Y_7[l]],label="",st=:scatter, ms=m_size, msw=0,mc=m_color, framestyle=:none,aspect_ratio=:equal)
        end
    end

    savefig(output_spin_rel_filename)
    return plt_c
end


#この関数を利用するときはファイル名ではなく、ディレクトリ名を入力することに注意
#ディレクトリ名には最後のスラッシュを含めないことに注意
function summary_plot_27site_kagome(marker_magnification, line_magnification, input_dir_name, output_dir_name, J_r, J_g, J_b, title_name)
    anim = Animation()
    
    input_MHdata_filename = input_dir_name * "/MHdata.csv"
    #サイト間のbondの強さをplotする
    plt_lattice = kagome_spinrel.plot_27site_kagome_lattice(J_r, J_g, J_b)
    for M_index=14:27
        #ディレクトリ名をもとに各データの入力ファイル名を指定する
        input_sz_filename = input_dir_name * "/sz/sz_" * string(M_index) * "_upstate.csv"
        input_szz_filename = input_dir_name * "/szz/szz_" * string(M_index) * "_nn_list.csv"
        input_sxx_filename = input_dir_name * "/sxx/sxx_" * string(M_index) * "_nn_list.csv"

        #ディレクトリ名をもとに各データの出力ファイル名を指定する
        output_spin_rel_filename = output_dir_name * "/spin_rel/spin_rel_" * string(M_index) * ".png"
        output_summary_filename = output_dir_name * "/summary/summary_" * string(M_index) * ".png"

        #plotの実行
        plt_c = kagome_spinrel.plot_27site_NN_spin_rel(input_szz_filename,input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_spin_rel_filename)
        plt_MHcurve = kagome_spinrel.plot_27site_MH_curve_with_Highrigt(input_MHdata_filename, M_index)
        
        
        M=[L"1/27",L"1/9",L"5/27",L"7/27",L"1/3 ",L"11/27",L"13/27",L"5/9 ",L"17/27",L"19/27",L"7/9 ",L"23/27 ",L"25/27",L"1  ",L"1  "]
        title=plot(title=title_name * M[M_index - 13],grid=false, titleposition =:left,showaxis=false,titleframe=:box)
        plot(title=title_name * M[M_index - 13],grid=false, titleposition =:left,showaxis=false,titleframe=:box)
        l = @layout[a{0.01h}; b c; d]
        plt_sum = plot(title ,plt_MHcurve, plt_lattice,plt_c, layout=l,size=(1640,1640),margin=10mm, plot_titlevspan=0.01)
        savefig(output_summary_filename)
        frame(anim, plt_sum)        
    end
    output_animation_filename = output_dir_name * "/animation.gif"
    gif(anim, output_animation_filename, fps = 1)
end

    function summary_plot_36site_kagome(M_start,M_end,marker_magnification, line_magnification, input_dir_name, output_dir_name, J_r, J_g, J_b, title_name)
        anim = Animation()
    
        input_MHdata_filename = input_dir_name * "/MHdata.csv"
        plt_lattice = kagome_spinrel.plot_36site_kagome_lattice(J_r, J_g, J_b)
        for M_index=M_start:M_end
            #ディレクトリ名をもとに各データの入力ファイル名を指定する
            input_sz_filename = input_dir_name * "/sz/sz_" * string(M_index) * "_upstate.csv"
            input_szz_filename = input_dir_name * "/szz/szz_" * string(M_index) * "_nn_list.csv"
            input_sxx_filename = input_dir_name * "/sxx/sxx_" * string(M_index) * "_nn_list.csv"

            #ディレクトリ名をもとに各データの出力ファイル名を指定する
            output_spin_rel_filename = output_dir_name * "/spin_rel/spin_rel_" * string(M_index) * ".png"
            output_summary_filename = output_dir_name * "/summary/summary_" *string(M_index) *".png"

            #plotの実行
            plt_c = kagome_spinrel.plot_36site_NN_spin_rel(input_szz_filename,input_sxx_filename,input_sz_filename ,marker_magnification,line_magnification,output_spin_rel_filename)
            plt_MHcurve = kagome_spinrel.plot_36site_MH_curve_with_Highrigt(input_MHdata_filename, M_index)
        
        
            M=[L"0",L"1/18",L"1/9",L"1/6",L"2/9",L"5/18",L"1/3 ",L"7/18",L"4/9",L"1/2",L"5/9",L"11/18",L"2/3",L"13/18",L"7/9 ",L"15/18",L"8/9 ",L"17/18",L"1  ",L"1"]
            plt_title=plot(title=title_name * M[M_index - 17], grid=false, titleposition =:left,showaxis=false)
            l = @layout[a{0.01h}; b c;d]
            plt_sum = plot(plt_title ,plt_MHcurve, plt_lattice,plt_c, layout=l,size=(1640,1640))
            savefig(output_summary_filename)
            frame(anim, plt_sum)        
        end
        output_animation_filename = output_dir_name * "/summary_animation.gif"
        gif(anim, output_animation_filename, fps = 1)
    end

export plot_36site_kagome_lattice
export plot_27site_MH_curve_with_Highrigt
export plot_36site_MH_curve_with_Highrigt
export plot_27site_szz_NNcorellation
export plot_27site_sxx_NNcorellation
export plot_36site_szz_NNcorellation
export plot_36site_sxx_NNcorellation
export summary_plot_27site_kagome
export summary_plot_36site_kagome
end
