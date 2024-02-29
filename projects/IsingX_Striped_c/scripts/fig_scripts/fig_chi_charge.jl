
########################
## charge susceptibility
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{charge}} $", ylabelsize=30,
    xlabelsize=30, title= my_title
)


df_LT=sort(unique(df[:,[:T,:L]]),:T, rev=true)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("compress")] , 
        dfs[!,Symbol("Δcompress")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("compress")] ;  
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = L"$β=%$(1/Ts)$", markersize=10
    )

    CairoMakie.scatterlines!(top, rangeU, χ_charge_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
        
end

axislegend( position=(1, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_charge.png"), fig)
end


########################
## charge-XX susceptibility
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{charge}}^{xx} $", ylabelsize=30,
    xlabelsize=30, title= my_title
)


df_LT=sort(unique(df[:,[:T,:L]]),:T, rev=true)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("CDSxx_00")] , 
        dfs[!,Symbol("ΔCDSxx_00")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("CDSxx_00")] ;  
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = L"$β=%$(1/Ts)$", markersize=10
    )

    CairoMakie.scatterlines!(top, rangeU, χ_charge_XX_analytic.(rangeU, 1/Ts, μ0_plot, Ls) *Ls^2;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
        
end

axislegend( position=(1, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_charge_xx.png"), fig)
end



################
## B1 charge density wave, (nematic proxy) figure
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{charge}}^{B_1} [\mathbf{Q}=(0,0)]$", 
    title=my_title, xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("B1_CDS_00")] , 
        dfs[!,Symbol("ΔB1_CDS_00")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("B1_CDS_00")]  ;  
        marker = :xcross,  color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))",  markersize=10
    )
    CairoMakie.scatterlines!(top, rangeU, χ_charge_B1_analytic.(rangeU, 1/Ts, μ0_plot, Ls) ./Ts;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )

end
axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_charge_B1_00.png"), fig)




################
## A1' charge density wave, (A1' bilinear proxy) figure
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{charge}}^{A_1^{\prime}} [\mathbf{Q}=(π,π)]$", 
    title=my_title, xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0,1))
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("CDS_ππ")] , 
        dfs[!,Symbol("ΔCDS_ππ")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("CDS_ππ")]  ;  
        marker = :xcross,  color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))", markersize=10
    )
    CairoMakie.scatterlines!(top, rangeU, χ_charge_A1P_analytic.(rangeU, 1/Ts, μ0_plot, Ls) *Ls^2;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_charge_A1p_ππ.png"), fig)



