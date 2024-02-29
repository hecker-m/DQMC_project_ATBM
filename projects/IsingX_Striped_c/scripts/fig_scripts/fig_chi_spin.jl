

########################
## spin susceptibility, Q=(0,0)
########################
q0=["0", "0"]
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, title=my_title
)


for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("SDS_Mx_x" * "_$(q0[1])$(q0[2])")], 
            dfs[!,Symbol("ΔSDS_Mx_x" * "_$(q0[1])$(q0[2])")]; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDS_Mx_x" * "_$(q0[1])$(q0[2])")] ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=20)
end
axislegend( position=(0.04,0.97))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q($(q0[1]),$(q0[2])).png"), fig)

########################
## spin susceptibility, Q=(π,0) + Q=(0,π)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{A_1}$", ylabelsize=30,
    xlabelsize=30, title=my_title2
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("SDS_A1_Mx_x")] , 
        dfs[!,Symbol("ΔSDS_A1_Mx_x")] ,
        df_ϕs[!,Symbol("ΔSDS_A1_00")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDS_A1_Mx_x")] ;
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")] , 
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax], 
         markersize=10
    )      

    CairoMakie.scatterlines!(top, rangeU, χspin_A1_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    
end

axislegend( position=(0.04,0.97))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_A1.png"), fig)

########################
## spin susceptibility, Q=(π,0) + Q=(0,π), as function of temperature
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$χ_{\mathrm{spin}}^{A_1}$", ylabelsize=30,
    xlabelsize=30, title=my_title2
)


for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  


    eb=errorbars!(top, dfs[!,:T], dfs[!,Symbol("SDS_A1_Mx_x")] , dfs[!,Symbol("ΔSDS_A1_Mx_x")] ,
        df_ϕs[!,Symbol("ΔSDS_A1_00")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("SDS_A1_Mx_x")]  ;
        marker = :xcross, color = colorschemes[:tab20][1-idx/num_T_points], 
        label = "U=$(Us)", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("SDS_A1_00")] , 
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], 
        markersize=10
    )         

    CairoMakie.scatterlines!(top, rangeT, χspin_A1_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
    marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
    color = colorschemes[:tab20][1-idx/num_T_points]
    )
end
axislegend(top, position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_A1_T.png"), fig)


################
## spin correlation function, Q=(π,0) + Q=(0,π)
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$S_{\mathrm{spin}}^{A_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=my_title2
)


for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("SDC_A1_Mx_x")] , 
        dfs[!,Symbol("ΔSDC_A1_Mx_x")] , 
        df_ϕs[!,Symbol("ΔSDC_A1_00")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDC_A1_Mx_x")]   ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("SDC_A1_00")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=10
    ) 
    CairoMakie.scatterlines!(top, rangeU, Sspin_A1_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
    marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
    color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )

end
axislegend( position=(0, 1))

################
## spin correlation function S(2), Q=(π,0) + Q=(0,π)
################
mid=Axis(fig[2, 1],  ylabel=L"$S_{\mathrm{spin}}^{(2), A_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(mid, dfs[!,:U], df_ϕs[!,Symbol("Sspin_2_A1_xx")] , 
        df_ϕs[!,Symbol("ΔSspin_2_A1_xx")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("Sspin_2_A1_xx")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=10
    ) 
    CairoMakie.scatterlines!(mid, rangeU, Sspin2_A1_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )

end

########################
## Magnetic Order parameter
########################
bottom=Axis(fig[3, 1], xlabel=L"$U/t$", ylabel=L"$\langle \left| \left( \mathbf{M}_{\mathbf{Q}_1}, 
    \mathbf{M}_{\mathbf{Q}_2}   \right)  \right| \rangle $", ylabelsize=30,
    xlabelsize=30, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_OP))  
 
    eb=errorbars!(bottom, df_OPs[!,:U], df_OPs[!,Symbol("Mx_X_OP34")], df_OPs[!,Symbol("ΔMx_X_OP34")],
        df_ϕ_OPs[!,Symbol("ΔMx_OP")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_OPs[!,:U], df_OPs[!,Symbol("Mx_X_OP34")] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10
    ) 
        
    # CairoMakie.scatter!(top, df_ϕ_OPs[!,:U], df_ϕ_OPs[!,Symbol("Mx_OP")]  ;  
    #     marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15
    # )         
end

display(fig)
CairoMakie.save(joinpath(p, "corr_spin_A1.png"), fig)



################
## spin correlation function, Q=(π,0) + Q=(0,π), as a function of T
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$S_{\mathrm{spin}}^{A_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=my_title2
)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(top, dfs[!,:T], dfs[!,Symbol("SDC_A1_Mx_x")] , 
        dfs[!,Symbol("ΔSDC_A1_Mx_x")] , 
        df_ϕs[!,Symbol("ΔSDC_A1_00")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("SDC_A1_Mx_x")]  ;  marker = :xcross,
        color = color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10
    )

    CairoMakie.scatter!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("SDC_A1_00")] ;  
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=10
    ) 
    CairoMakie.scatterlines!(top, rangeT, Sspin_A1_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
    marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
    color = colorschemes[:tab20][1-idx/num_T_points]
    )
end
axislegend( position=(1, 1))

################
## spin correlation function S(2), Q=(π,0) + Q=(0,π), as a function of T
################
mid=Axis(fig[2, 1],  ylabel=L"$S_{\mathrm{spin}}^{(2), A_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(mid, dfs[!,:T], df_ϕs[!,Symbol("Sspin_2_A1_xx")] , 
        df_ϕs[!,Symbol("ΔSspin_2_A1_xx")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:T], df_ϕs[!,Symbol("Sspin_2_A1_xx")] ;  
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=10
    ) 
    CairoMakie.scatterlines!(mid, rangeT, Sspin2_A1_XX_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][1-idx/num_T_points]
    )
end


########################
## Magnetic Order parameter, as a function of temperature
########################
bottom=Axis(fig[3, 1], xlabel=L"$T/t$", ylabel=L"$\langle \left| \left( \mathbf{M}_{\mathbf{Q}_1}, 
\mathbf{M}_{\mathbf{Q}_2}   \right)  \right| \rangle $", xlabelsize=30, ylabelsize=30,
xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

mkey="Mx_X_OP34"

df_OP_LU=sort(unique(df_OP[:,[:U,:L]]),:U)
for (idx, (Us, Ls)) in enumerate(zip(df_OP_LU[:,:U],df_OP_LU[:,:L]))
    num_T_points=length(df_OP_LU[:,:U]);
    df_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_OP))  
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ_OP))  


    eb=errorbars!(bottom, df_OPs[!,:T], df_OPs[!,Symbol(mkey)], df_OPs[!,Symbol("Δ" * mkey)]; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_OPs[!,:T], df_OPs[!,Symbol(mkey)] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], 
        label = "U=$(round(Us, digits=2)), L=$Ls", markersize=10
    )  
        
    # CairoMakie.scatter!(bottom, df_ϕ_OPs[!,:T], df_ϕ_OPs[!,Symbol("Mx_OP")]  ;  
    #     marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=15
    # )      
end
#axislegend( position=(1,1))

display(fig)
CairoMakie.save(joinpath(p, "corr_spin_A1_T.png"), fig)





