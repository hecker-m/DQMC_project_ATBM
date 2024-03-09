

################
## B1' bilinear susceptibility
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{B_1^{\prime}}$", 
    title=my_title, xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)


for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("B1p_dQ_S_00")] , 
        dfs[!,Symbol("ΔB1p_dQ_S_00")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("B1p_dQ_S_00")]  ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = L"$β=%$(1/Ts)$", markersize=10
    )
end

axislegend( position=(1, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_bilinear_B1p_00.png"), fig)
end


################
## B1' bilinear correlation function
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1],  xlabel= L"$U/t$", ylabel=L"$S_{\mathrm{bilinear}}^{B_1^{\prime}}$", 
    title=my_title, xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("B1p_dQ_C_00")] , 
        dfs[!,Symbol("ΔB1p_dQ_C_00")]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("B1p_dQ_C_00")]  ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = L"$β=%$(1/Ts)$", markersize=10
    )
end
axislegend( position=(0, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "corr_bilinear_B1p_00.png"), fig)
end
