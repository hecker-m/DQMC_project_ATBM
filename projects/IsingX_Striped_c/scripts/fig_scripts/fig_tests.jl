

################
## Nematic susceptibility vs. magnetic susceptibility for various U values
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,\pi)] +
    χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,0)]$", 
    ylabel=L"$χ_{\mathrm{nematic}}^{B_1} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30,     title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  


    eb=errorbars!(df_ϕs[!,Symbol("SDS_A1_00")], df_ϕs[!,Symbol("NemS_00")],
        df_ϕs[!,Symbol("ΔNemS_00")]; direction = :y, linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)

    eb_x=errorbars!(df_ϕs[!,Symbol("SDS_A1_00")], df_ϕs[!,Symbol("NemS_00")],
        df_ϕs[!,Symbol("ΔSDS_A1_00")] ; direction = :x, linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb_x, 0, 0, -0.5)

    CairoMakie.scatterlines!(top, df_ϕs[!,Symbol("SDS_A1_00")] , 
        df_ϕs[!,Symbol("NemS_00")],
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], 
         markersize=15, label = "U=$(Us)")         
end
axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_nematic_vs_chi_spin.png"), fig)



################
## Zoom: Nematic susceptibility vs. magnetic susceptibility for various U values
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,\pi)] +
    χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,0)]$", 
    ylabel=L"$χ_{\mathrm{nematic}}^{B_1} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30,     title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle", 
    limits=((0,30), (0,80)))

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    if Us >0.8
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  


    eb=errorbars!(df_ϕs[!,Symbol("SDS_A1_00")], df_ϕs[!,Symbol("NemS_00")],
        df_ϕs[!,Symbol("ΔNemS_00")]; direction = :y, linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)

    eb_x=errorbars!(df_ϕs[!,Symbol("SDS_A1_00")], df_ϕs[!,Symbol("NemS_00")],
        df_ϕs[!,Symbol("ΔSDS_A1_00")] ; direction = :x, linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb_x, 0, 0, -0.5)

    CairoMakie.scatterlines!(top, df_ϕs[!,Symbol("SDS_A1_00")] , 
        df_ϕs[!,Symbol("NemS_00")],
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], 
         markersize=15, label = "U=$(Us)")    
    end     
end
axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_nematic_vs_chi_spin_zoom.png"), fig)
