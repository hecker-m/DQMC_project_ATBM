########################
## phase stiffness
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$\rho_{s}\;\cdot  β\,\frac{\pi}{2} $", ylabelsize=30,
    xlabelsize=30, title=my_title
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
 
    eb=errorbars!(top, dfs[!,:U], dfs[!,:ρ_s] *π/(2Ts), dfs[!,:Δρ_s] *π/(2Ts); linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, 0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:ρ_s] *π/(2Ts);  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = L"$β=%$(1/Ts)$", markersize=10
    )
end
axislegend( position=(1, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "phase_stiffness.png"), fig)
end

########################
## check ΛL ≈ kx
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$\Lambda_{L} - \langle k_x \rangle $", ylabelsize=30,
    xlabelsize=30, title=my_title
)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
 
    eb=errorbars!(top, dfs[!,:U], dfs[!,:ΛL] .+ dfs[!,:k_x] , dfs[!,:ΔΛL] .+ dfs[!,:Δk_x]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax]
    )
    CairoMakie.translate!(eb, 0, 0, 0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:ΛL] .+ dfs[!,:k_x];  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = L"$β=%$(1/Ts)$", markersize=10
    )
end
#axislegend( position=(0, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "phase_stiffness_kx_check.png"), fig)
end
