

################
## B1' bilinear susceptibility
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{B_1^{\prime}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("B1p_dQ_S_00")] , 
            dfs[!,Symbol("ΔB1p_dQ_S_00")], df_ϕs[!,Symbol("ΔB1p_dQ_S_00")] /3; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("B1p_dQ_S_00")]  ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=ms)

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("B1p_dQ_S_00")] /3;  
            marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=ms_ϕ) 
end

axislegend( position=(1, 1))
display(fig)
CairoMakie.save(joinpath(p, "chi_bilinear_B1p_00.png"), fig)



################
## B1' bilinear correlation function
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1],  xlabel= L"$U/t$", ylabel=L"$\mathrm{corr:\,}S_{\mathrm{bilinear}}^{B_1^{\prime}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
 
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("B1p_dQ_C_00")] , 
            dfs[!,Symbol("ΔB1p_dQ_C_00")], df_ϕs[!,Symbol("ΔB1p_dQ_C_00")] /3; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("B1p_dQ_C_00")]  ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=ms)

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("B1p_dQ_C_00")] /3;  
            marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=ms_ϕ) 
    
end
axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "corr_bilinear_B1p_00.png"), fig)
