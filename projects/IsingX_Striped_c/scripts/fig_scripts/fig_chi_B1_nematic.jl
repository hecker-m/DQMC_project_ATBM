

################
## B1 nematic susceptibility
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{B_1}$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, 
    limits= (nothing, nothing), title=my_title2
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
    
    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("NemS_X")] , 
        dfs[!,Symbol("ΔNemS_X")], df_ϕs[!,Symbol("ΔNemS")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("NemS_X")]  ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = L"$β=%$(1/Ts)$", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("NemS")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=10
    )         
    
    CairoMakie.scatterlines!(top, rangeU, χbil_B1_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )


end

axislegend( position=(0, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_nematic_B1_00.png"), fig)
end

################
## B1 nematic susceptibility as a function of T
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{B_1} $", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, 
    limits= (nothing, nothing),title=my_title2
)


for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("NemS_X")] , 
        dfs[!,Symbol("ΔNemS_X")], df_ϕs[!,Symbol("ΔNemS")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("NemS_X")]  ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("NemS")];  
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=10
    )            
    
    CairoMakie.scatterlines!(top, rangeT, χbil_B1_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][1-idx/num_T_points]
    )
end

axislegend( position=(1, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_nematic_B1_00_T.png"), fig)
end

################
## B1 nematic correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$S_{\mathrm{bilinear}}^{B_1} $", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=my_title2
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("NemC_X")] , 
        dfs[!,Symbol("ΔNemC_X")], df_ϕs[!,Symbol("ΔNemC")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("NemC_X")]  ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = L"$β=%$(1/Ts)$", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("NemC")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=10
    ) 
    CairoMakie.scatterlines!(top, rangeU, Sbil_B1_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
end
axislegend(top, position=(0, 1))

################
## B1 nematic correlation function S(2)
################
mid=Axis(fig[2, 1],  ylabel=L"$S_{\mathrm{bilinear}}^{(2), B_1} $", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(mid, dfs[!,:U], df_ϕs[!,Symbol("Sbil_2_B1_xx")] , 
        df_ϕs[!,Symbol("ΔSbil_2_B1_xx")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_B1_xx")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=10
    ) 
    CairoMakie.scatterlines!(mid, rangeU, Sbil2_B1_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
end

################
## B1 nematic order parameter 
################
bottom=Axis(fig[3, 1], xlabel=L"$U/t$", ylabel=L"$\left| \Phi^{B_1} \right|$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_OP))  
 
    eb=errorbars!(bottom, df_OPs[!,:U], df_OPs[!,Symbol("B1_OP")], df_OPs[!,Symbol("ΔB1_OP")];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_OPs[!,:U], df_OPs[!,Symbol("B1_OP")] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = L"$β=%$(1/Ts)$", markersize=10
    ) 
        
    # CairoMakie.scatter!(top, df_ϕ_OPs[!,:U], df_ϕ_OPs[!,Symbol("B1_OP")] ;  
    #     marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15
    # ) 
end
#axislegend( position=(0, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "corr_nematic_B1_00.png"), fig)
end

################
## B1 nematic correlation function as function of temperature
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], ylabel=L"$S_{\mathrm{bilinear}}^{B_1}$", 
    xlabelsize=30, ylabelsize=30,xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=my_title2
)


for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  
 
eb=errorbars!(top, dfs[!,:T], dfs[!,Symbol("NemC_X")] , 
        dfs[!,Symbol("ΔNemC_X")], df_ϕs[!,Symbol("ΔNemC")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("NemC_X")]  ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("NemC")] , marker = '□',
        color = colorschemes[:tab20][1-idx/num_T_points], markersize=10
    )         
    
    CairoMakie.scatterlines!(top, rangeT, Sbil_B1_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][1-idx/num_T_points]
    )
end
axislegend( position=(1, 1))

################
## B1 nematic correlation function S(2) as function of temperature
################
mid=Axis(fig[2, 1], ylabel=L"$S_{\mathrm{bilinear}}^{(2), B_1}$", 
    xlabelsize=30, ylabelsize=30,xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)


for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  
 
    eb=errorbars!(mid, dfs[!,:T], df_ϕs[!,Symbol("Sbil_2_B1_xx")] , 
        df_ϕs[!,Symbol("ΔSbil_2_B1_xx")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)

    CairoMakie.scatterlines!(mid, df_ϕs[!,:T], df_ϕs[!,Symbol("Sbil_2_B1_xx")] , marker = '□',
        color = colorschemes[:tab20][1-idx/num_T_points], markersize=10
    )         
    
    CairoMakie.scatterlines!(mid, rangeT, Sbil2_B1_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][1-idx/num_T_points]
    )
end

################
## B1 nematic order parameter as function of temperature
################
bottom=Axis(fig[3, 1], xlabel=L"$T/t$", ylabel=L"$\left| \Phi^{B_1} \right|$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Us, Ls)) in enumerate(zip(df_OP_LU[:,:U],df_OP_LU[:,:L]))
    num_T_points=length(df_OP_LU[:,:U]);
    df_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ_OP))  
 
    eb=errorbars!(bottom, df_OPs[!,:T], df_OPs[!,Symbol("B1_OP")], 
    df_OPs[!,Symbol("ΔB1_OP")], df_ϕ_OPs[!,Symbol("ΔB1_OP")];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_OPs[!,:T], df_OPs[!,Symbol("B1_OP")] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10
    ) 
        
    # CairoMakie.scatter!(top, df_ϕ_OPs[!,:T], df_ϕ_OPs[!,Symbol("B1_OP")] ;  
    #     marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=15
    # ) 
end
#axislegend( position=(0, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "corr_nematic_B1_00_T.png"), fig)
end

