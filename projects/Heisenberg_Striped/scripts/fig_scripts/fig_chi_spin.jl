

########################
## spin susceptibility, Q=(0,0)
########################
q0=["0", "0"]
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("SDS_Mx_x" * "_$(q0[1])$(q0[2])")], 
            dfs[!,Symbol("ΔSDS_Mx_x" * "_$(q0[1])$(q0[2])")]; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDS_Mx_x" * "_$(q0[1])$(q0[2])")] ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=ms)
end
axislegend( position=(0.04,0.97))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q($(q0[1]),$(q0[2])).png"), fig)

########################
## spin susceptibility, Q=(π,0) + Q=(0,π)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,\pi)] +
    χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,0)]$", ylabelsize=30,
    xlabelsize=30, title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(dfs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")]*Ls^2 /3, 
            dfs[!,Symbol("ΔSDS_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_x_π0")],
            df_ϕs[!,Symbol("ΔSDS_A1_00")]*Ls^2 /3; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")] ;
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))", markersize=ms)

    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")]*Ls^2 /3, 
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax], 
         markersize=ms_ϕ)      
end

axislegend( position=(0.04,0.97))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q_0π_Plus_π0.png"), fig)

########################
## spin susceptibility, Q=(π,0) + Q=(0,π), as function of temperature
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,\pi)] +
    χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,0)]$", ylabelsize=30,
    xlabelsize=30, title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")


CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  


    eb=errorbars!(dfs[!,:T], df_ϕs[!,Symbol("SDS_A1_00")]*Ls^2 /3, 
            dfs[!,Symbol("ΔSDS_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_x_π0")],
            df_ϕs[!,Symbol("ΔSDS_A1_00")]*Ls^2 /3; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:T], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")] ;
        marker = :xcross, color = colorschemes[:tab20][1-idx/num_T_points], 
        label = "U=$(Us)", markersize=ms)

    CairoMakie.scatterlines!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("SDS_A1_00")]*Ls^2 /3 , 
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], 
         markersize=ms_ϕ)         
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q_0π_Plus_π0_T.png"), fig)


########################
## spin susceptibility, Q=(π,0) - Q=(0,π)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,\pi)] -
    χ_{\mathrm{spin}}^{} [\mathbf{Q}=(\pi,0)]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=8")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("SDS_Mx_x_0π")] .- dfs[!,Symbol("SDS_Mx_x_π0")], 
            dfs[!,Symbol("ΔSDS_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_x_π0")]; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDS_Mx_x_0π")] .- dfs[!,Symbol("SDS_Mx_x_π0")] ;
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))", markersize=ms)
end
axislegend( position=(0.04,0.97))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q_0π_Minus_π0.png"), fig)


########################
## Magnetic Order parameter
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$\langle \left| \left( \mathbf{M}_{\mathbf{Q}_1}, 
    \mathbf{M}_{\mathbf{Q}_2}   \right)  \right| \rangle $", ylabelsize=30,
    xlabelsize=30, title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), (0, 0, -0.8))

for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_OP))  
 
    eb=errorbars!(df_OPs[!,:U], df_OPs[!,Symbol("Mx_X_OP34")], df_OPs[!,Symbol("ΔMx_X_OP34")],
        df_ϕ_OPs[!,Symbol("ΔMx_OP")] /3; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:U], df_OPs[!,Symbol("Mx_X_OP34")] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=ms) 
        
    CairoMakie.scatter!(top, df_ϕ_OPs[!,:U], df_ϕ_OPs[!,Symbol("Mx_OP")] /3 ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=ms_ϕ)         
end
# CairoMakie.translate!(poly!(Point2f[(1.5, 0), (1.5, end), (2.4, end), (2.4, 0)], color = :gray, strokewidth = 0), (0, 0, -0.8))
axislegend( position=(1,1))

display(fig)
CairoMakie.save(joinpath(p, "OP_Magnetic.png"), fig)

## as a function of temperature
##
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel="T/t", ylabel="⟨|[M(Q₁), M(Q₂)]|⟩", xlabelsize=30, ylabelsize=30,
xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

mkey="Mx_X_OP34"

df_OP_LU=sort(unique(df_OP[:,[:U,:L]]),:U)
for (idx, (Us, Ls)) in enumerate(zip(df_OP_LU[:,:U],df_OP_LU[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_OP))  
    eb=errorbars!(df_OPs[!,:T], df_OPs[!,Symbol(mkey)], df_OPs[!,Symbol("Δ" * mkey)]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][idx])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:T], df_OPs[!,Symbol(mkey)] ;  marker = :xcross,
        color = colorschemes[:tab20][idx], label = "U=$(round(Us, digits=2)), L=$Ls", markersize=ms)       
end
axislegend( position=(1,1))

display(fig)
CairoMakie.save(joinpath(p, "OP_Magnetic_T.png"), fig)


