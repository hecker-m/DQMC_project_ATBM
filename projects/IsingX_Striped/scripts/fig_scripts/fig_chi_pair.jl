df_LT=sort(unique(df[:,[:T,:L]]),:T)

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
axislegend( position=(1,1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_00_00.png"), fig)

########################
## pairing susceptibility, 00 Q=(π0)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{00} [\mathbf{Q}=(π,0)]$", ylabelsize=30,
    xlabelsize=30, xticklabelsize=20, yticklabelsize=20, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_s_π0], dfs[!,:ΔPDS_s_π0]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:PDS_s_π0] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$(1/Ts)", markersize=20)        
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_00_π0.png"), fig)

########################
## pairing susceptibility,  ZZ
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{zz} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, xticklabelsize=20, yticklabelsize=20, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_spm_00], dfs[!,:ΔPDS_spm_00]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:PDS_spm_00] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$(1/Ts)", markersize=20)        
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_ZZ_00.png"), fig)

########################
## pairing susceptibility,  ZZ Q=(π0)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{zz} [\mathbf{Q}=(π,0)]$", ylabelsize=30,
    xlabelsize=30, xticklabelsize=20, yticklabelsize=20, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_spm_π0], dfs[!,:ΔPDS_spm_π0]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,:PDS_spm_π0] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$(1/Ts)", markersize=20)        
end

axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_ZZ_π0.png"), fig)

########################
## pairing susceptibility, states XX 
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{xx} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, xticklabelsize=20, yticklabelsize=20, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)


for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    nTs=length(df_LT[!, :T]);
    i1=1+(idx-1)*2
    i2=2+(idx-1)*2
    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_XX_00], dfs[!,:ΔPDS_XX_00]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:rainbow][i1/(2nTs)])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:PDS_XX_00] ;  marker = :xcross,
        color = colorschemes[:rainbow][i1/(2nTs)], label = L"\mathbf{Q}=(0,0), β=%$(1/Ts)", markersize=20)

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_XX_π0], dfs[!,:ΔPDS_XX_π0]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:rainbow][i2/(2nTs)])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:PDS_XX_π0] ;  marker = :xcross,
        color = colorschemes[:rainbow][i2/(2nTs)], label = L"\mathbf{Q}=(π,0), β=%$(1/Ts)", markersize=20)
       
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_XX.png"), fig)
########################
## pairing susceptibility, states YYzz, Q=00
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{YYzz} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, xticklabelsize=20, yticklabelsize=20, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    nTs=length(df_LT[!, :T]);
    i3=1+(idx-1)*1

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_YYzz_00], dfs[!,:ΔPDS_YYzz_00]; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:rainbow][i3/(1nTs)])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:PDS_YYzz_00] ;  marker = :xcross,
        color = colorschemes[:rainbow][i3/(1nTs)], label = "β=$(1/Ts)", markersize=20)     
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_YYzz_00.png"), fig)

########################
## pairing susceptibility, states YYzz Q=(π0)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{pair}}^{YYzz} [\mathbf{Q}=(π,0)]$", ylabelsize=30,
    xlabelsize=30, xticklabelsize=20, yticklabelsize=20, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    nTs=length(df_LT[!, :T]);
    i4=1+(idx-1)*1

    eb=errorbars!(dfs[!,:U], dfs[!,:PDS_YYzz_π0], dfs[!,:ΔPDS_YYzz_π0]; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:rainbow][i4/(1nTs)])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:PDS_YYzz_π0] ;  marker = :xcross,
        color = colorschemes[:rainbow][i4/(1nTs)], label = "β=$(1/Ts)", markersize=20)        
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_pair_YYzz_π0.png"), fig)


########################
## SC bilinear order parameter Δ
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel="U/t", ylabel="⟨ΔᶻʸΔᶻʸ⟩", xlabelsize=30, ylabelsize=30,
xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

mkey="Δ_Zy_bil_OP1"

df_OP_LT=sort(unique(df_OP[:,[:T,:L]]),:T)
for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP))  
    eb=errorbars!(df_OPs[!,:U], df_OPs[!,Symbol(mkey)], df_OPs[!,Symbol("Δ" * mkey)]; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:tab20][idx])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:U], df_OPs[!,Symbol(mkey)] ;  marker = :xcross,
        color = colorschemes[:tab20][idx], label = "T=$(round(Ts, digits=2)), L=$Ls", markersize=20)       
end
axislegend( position=(1,1))

display(fig)
CairoMakie.save(joinpath(p, "OP_SC_bil_ZZyy.png"), fig)


fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel="T/t", ylabel="⟨ΔᶻʸΔᶻʸ⟩", xlabelsize=30, ylabelsize=30,
xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

df_OP_LU=sort(unique(df_OP[:,[:U,:L]]),:U)
for (idx, (Us, Ls)) in enumerate(zip(df_OP_LU[:,:U],df_OP_LU[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_OP))  
    eb=errorbars!(df_OPs[!,:T], df_OPs[!,Symbol(mkey)], df_OPs[!,Symbol("Δ" * mkey)]; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:tab20][idx])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:T], df_OPs[!,Symbol(mkey)] ;  marker = :xcross,
        color = colorschemes[:tab20][idx], label = "U=$(round(Us, digits=2)), L=$Ls", markersize=20)       
end
axislegend( position=(1,1))

display(fig)
CairoMakie.save(joinpath(p, "OP_SC_bil_ZZyy_T.png"), fig)