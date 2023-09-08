
########################
## pairing susceptibility, 00
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{00} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_s_00], dfs[!,:ΔPDS_s_00]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:PDS_s_00] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$(1/Ts)", markersize=20)        
end
axislegend( position=(0.8,1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_00_00.png"), fig)



########################
## pairing susceptibility,  ZZ
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{zz} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, xticklabelsize=20, yticklabelsize=20, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

χc_ref = FileIO.load(joinpath(p, "ScreenShots/screenshot_chi_pair.png"))
ip = image!(top, 0..3.05, -1.38..41, χc_ref'[:, end:-1:1],transparency=true)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_spm_00], dfs[!,:ΔPDS_spm_00]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,:PDS_spm_00] ;  marker = '□',
        color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$(1/Ts)", markersize=15)        
end
axislegend( position=(0.8, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_ZZ_00.png"), fig)




