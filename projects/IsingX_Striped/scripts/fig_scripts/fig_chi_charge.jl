
########################
## charge susceptibility
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{charge}} $", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2,0.3, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

df_LT=sort(unique(df[:,[:T,:L]]),:T, rev=true)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
    # eb=errorbars!(dfs[!,:U], dfs[!,:CDS], dfs[!,:Δ_CDS]; linewidth=2.5, whiskerwidth=10,
    #     color = colorschemes[:tab20][idx])
    # CairoMakie.translate!(eb, 0, 0, -0.5) 
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:CDS_00] ;  marker = :xcross,
    color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=20)
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_charge.png"), fig)