########################
## phase stiffness
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$\rho_{s}\;\cdot  β\,\frac{\pi}{2} $", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

# CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.0, 1.0], color = :gray,linewidth=0.2), 0, 0, 0.3)
ρs_ref = FileIO.load(joinpath(p_screenshot, "ScreenShots/screenshot_phasestiffness.png"))
ip = image!(top, 0..3.05, -0.34..1.45, ρs_ref'[:, end:-1:1],transparency=true)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
 
    eb=errorbars!(dfs[!,:U], dfs[!,:rho_s_00_4] *π/(2Ts), dfs[!,:Δrho_s_00_4] *π/(2Ts); linewidth=er_lw, 
    whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, 0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,:rho_s_00_4] *π/(2Ts);  marker = '□',
    color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=20)
end
axislegend( position=(0.8, 1))

display(fig)
CairoMakie.save(joinpath(p, "phase_stiffness.png"), fig)


########################
## check ΛL ≈ kx
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$\Lambda_{L} - \langle k_x \rangle $", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

# CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
# CairoMakie.translate!(hlines!(top, [0.0, 1.0], color = :gray,linewidth=0.2), 0, 0, 0.3)


for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
 
    eb=errorbars!(dfs[!,:U], dfs[!,:rho_s_00_1] .+ dfs[!,:kx_00] , dfs[!,:Δrho_s_00_1] .+ dfs[!,:Δkx_00]; linewidth=er_lw, 
    whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, 0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,:rho_s_00_1] .+ dfs[!,:kx_00];  marker = :xcross,
    color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=20)
end
axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "phase_stiffness_kx_check.png"), fig)