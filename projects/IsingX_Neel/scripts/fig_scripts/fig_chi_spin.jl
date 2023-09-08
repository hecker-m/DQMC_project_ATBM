########################
## spin susceptibility
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,\pi)]$", ylabelsize=30,
    xlabelsize=30, 
    title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")


CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(dfs[!,:U], dfs[!,:SDS_Mx_x_ππ] .* 4 .* dfs[!,:U].^2, 
        dfs[!,:ΔSDS_Mx_x_ππ] .* 4 .* dfs[!,:U].^2, 
        df_ϕs[!,:ΔSDS_Mx_x_ππ] .* 4 .* df_ϕs[!,:U].^2 .*Ls^2; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:SDS_Mx_x_ππ] .* 4 .* dfs[!,:U].^2;  
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$(1/Ts)", markersize=10) 
        
    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,:SDS_Mx_x_ππ].* 4 .* df_ϕs[!,:U].^2 .*Ls^2;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax], 
         markersize=15) 
end



axislegend( position=(0,1))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_ππ.png"), fig)


########################
## spin susceptibility, rescaled with L^(-7/4), β=10
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$L^{-7/4}\cdot χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,\pi)]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

χc_ref = FileIO.load(joinpath(p, "ScreenShots/screenshot_chi_spin.png"))
ip = image!(top, 0.5..2.0, -6..122, χc_ref'[:, end:-1:1],transparency=true)
CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    if Ts==0.1
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    eb=errorbars!(dfs[!,:U], dfs[!,:SDS_Mx_x_ππ] .* 4 .*Ls^(-7/4) .* dfs[!,:U].^2, 
        dfs[!,:ΔSDS_Mx_x_ππ] .* 4 .*Ls^(-7/4) .* dfs[!,:U].^2; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][0.02])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:SDS_Mx_x_ππ] .* 4 .*Ls^(-7/4) .* dfs[!,:U].^2;  
        marker = '□', color = colorschemes[:coolwarm][0.02], 
        label = "β=$(1/Ts), L=8", markersize=15)  
    end      
end
axislegend( position=(0,1))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_ππ_rescaled_b10.png"), fig)

########################
## spin susceptibility, rescaled with L^(-7/4), β=20
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$L^{-7/4}\cdot χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,\pi)]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

χc_ref = FileIO.load(joinpath(p, "ScreenShots/screenshot_chi_spin_b20.png"))
ip = image!(top, 0.5..2.0, -8..170, χc_ref'[:, end:-1:1],transparency=true)
CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    if Ts==0.05
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  

    eb=errorbars!(dfs[!,:U], dfs[!,:SDS_Mx_x_ππ] .* 4 .*Ls^(-7/4) .* dfs[!,:U].^2, 
        dfs[!,:ΔSDS_Mx_x_ππ] .* 4 .*Ls^(-7/4) .* dfs[!,:U].^2; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][0.02])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,:SDS_Mx_x_ππ] .* 4 .*Ls^(-7/4) .* dfs[!,:U].^2;  
        marker = '□', color = colorschemes[:coolwarm][0.02], 
        label = "β=$(1/Ts), L=8", markersize=15)  
    end      
end
axislegend( position=(0,1))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_ππ_rescaled_b20.png"), fig)