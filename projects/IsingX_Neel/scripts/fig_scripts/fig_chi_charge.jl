########################
## charge susceptibility
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{charge}} $", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2,0.3, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)
χc_ref = FileIO.load(joinpath(p, "ScreenShots/screenshot_chi_c2.png"))
ip = image!(top, 0..3.05, -0.03..0.55, χc_ref'[:, end:-1:1],transparency=true)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
 
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,:CDS_00] ;  marker = '□',
    color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=20)
end
axislegend( position=(0.8, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_charge.png"), fig)