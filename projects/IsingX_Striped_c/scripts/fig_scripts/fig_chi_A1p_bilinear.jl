
################
## A1' bilinear susceptibility
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{A_1^{\prime}}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsize=20, yticklabelsize=20, 
    limits= (nothing, nothing), title=my_title2
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
 
    eb=errorbars!(dfs[!,:U], df_ϕs[!,Symbol("A1p_dQ_S")], 
        dfs[!,Symbol("ΔA1p_dQ_S_X")], df_ϕs[!,Symbol("ΔA1p_dQ_S")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,Symbol("A1p_dQ_S_X")]  ;  
        marker = :xcross,  color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))", markersize=10
    )

    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_S")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax], markersize=10
    )         
 
    # CairoMakie.scatterlines!(top, rangeU, χbil_A1P_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
    #     marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
    #     color = colorschemes[:coolwarm][1-1/Ts/βmax]
    # )
end

axislegend( position=(0, 1))
display(fig)
CairoMakie.save(joinpath(p, "chi_bilinear_A1p_00.png"), fig)

################
## A1' bilinear susceptibility as a function of T
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{A_1^{\prime}} $", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=my_title2
)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(dfs[!,:T], df_ϕs[!,Symbol("A1p_dQ_S")]  , 
        dfs[!,Symbol("ΔA1p_dQ_S_X")], df_ϕs[!,Symbol("ΔA1p_dQ_S")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:T], dfs[!,Symbol("A1p_dQ_S_X")]  ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", 
        markersize=10
    )

    CairoMakie.scatterlines!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("A1p_dQ_S")] ;  
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=10
    )            
        
    # CairoMakie.scatterlines!(top, rangeT, χbil_A1P_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
    #     marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
    #     color = colorschemes[:tab20][1-idx/num_T_points]
    # )
end

axislegend( position=(1, 1))
display(fig)
CairoMakie.save(joinpath(p, "chi_bilinear_A1p_00_T.png"), fig)


################
## A1' bilinear correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$S_{\mathrm{bilinear}}^{A_1^{\prime}} $", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=my_title2
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
 
    eb=errorbars!(top, dfs[!,:U], df_ϕs[!,Symbol("A1p_dQ_C")] , 
        dfs[!,Symbol("ΔA1p_dQ_C_X")], df_ϕs[!,Symbol("ΔA1p_dQ_C")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,Symbol("A1p_dQ_C_X")]  ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10
    )

    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_C")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  
        markersize=10
    ) 
    CairoMakie.scatterlines!(top, rangeU, Sbil_A1P_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
end
axislegend( position=(0, 1))

################
## A1' bilinear correlation function S(2)
################
mid=Axis(fig[2, 1],  ylabel=L"$S_{\mathrm{bilinear}}^{(2),A_1^{\prime}} $", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
 
    eb=errorbars!(mid, dfs[!,:U], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] , 
        df_ϕs[!,Symbol("ΔSbil_2_A1p_xx")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  
        markersize=10
    ) 
    CairoMakie.scatterlines!(mid, rangeU, Sbil2_A1P_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
end

################
## A1' bilinear order parameter 
################
bottom=Axis(fig[3, 1], xlabel=L"$U/t$", ylabel=L"$\left| \Phi^{A_1^{\prime}} \right|$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_OP))  
 
    eb=errorbars!(bottom,df_OPs[!,:U], df_OPs[!,Symbol("A1p_OP")], 
    df_OPs[!,Symbol("ΔA1p_OP")], df_ϕ_OPs[!,Symbol("ΔA1p_OP")];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_OPs[!,:U], df_OPs[!,Symbol("A1p_OP")] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", 
        markersize=10
    ) 
        
    # CairoMakie.scatter!(bottom, df_ϕ_OPs[!,:U], df_ϕ_OPs[!,Symbol("A1p_OP")] ;  
    #     marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15
    # )
end
#axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "corr_bilinear_A1p_00.png"), fig)

################
## A1' bilinear correlation function as function of temperature
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], ylabel=L"$S_{\mathrm{bilinear}}^{A_1^{\prime}} $", 
    xlabelsize=30, ylabelsize=30,xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=my_title2
)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(top, dfs[!,:T], df_ϕs[!,Symbol("A1p_dQ_C")] , 
        dfs[!,Symbol("ΔA1p_dQ_C_X")], df_ϕs[!,Symbol("ΔA1p_dQ_C")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:T], dfs[!,Symbol("A1p_dQ_C_X")]  ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", 
        markersize=10
    )
    CairoMakie.scatterlines!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("A1p_dQ_C")] , 
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], markersize=10
    )  

    CairoMakie.scatterlines!(top, rangeT, Sbil_A1P_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][1-idx/num_T_points]
    )
end
axislegend( position=(1, 1))

################
## A1' bilinear correlation function as function of temperature
################
mid=Axis(fig[2, 1], ylabel=L"$S_{\mathrm{bilinear}}^{(2),A_1^{\prime}} $", 
    xlabelsize=30, ylabelsize=30,xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(mid, dfs[!,:T], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] , 
        df_ϕs[!,Symbol("ΔSbil_2_A1p_xx")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:T], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] , 
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], markersize=10
    )  

    CairoMakie.scatterlines!(mid, rangeT, Sbil2_A1P_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][1-idx/num_T_points]
    )
end

################
## A1' bilinear order parameter as function of temperature
################
bottom=Axis(fig[3, 1], xlabel=L"$T/t$", ylabel=L"$\left| \Phi^{A_1^{\prime}} \right|$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Us, Ls)) in enumerate(zip(df_OP_LU[:,:U],df_OP_LU[:,:L]))
    num_T_points=length(df_OP_LU[:,:U]);
    df_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ_OP))  
 
    eb=errorbars!(bottom, df_OPs[!,:T], df_OPs[!,Symbol("A1p_OP")], df_OPs[!,Symbol("ΔA1p_OP")], 
        df_ϕ_OPs[!,Symbol("ΔA1p_OP")]; linewidth=er_lw, whiskerwidth=10,
        color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_OPs[!,:T], df_OPs[!,Symbol("A1p_OP")] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10
    )
    # CairoMakie.scatter!(bottom, df_ϕ_OPs[!,:T], df_ϕ_OPs[!,Symbol("A1p_OP")] ;  
    #     marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=15
    # ) 
       
end
#axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "corr_bilinear_A1p_00_T.png"), fig)