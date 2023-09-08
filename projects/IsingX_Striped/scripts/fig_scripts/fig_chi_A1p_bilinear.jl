df_LT=sort(unique(df[:,[:T,:L]]),:T, rev=true)
df_LU=sort(unique(df[:,[:U,:L]]),:U)
df_OP_LT=sort(unique(df_OP[:,[:T,:L]]),:T)
df_OP_LU=sort(unique(df_OP[:,[:U,:L]]),:U)

################
## A1' charge density wave, (A1' bilinear proxy) figure
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{charge}}^{A_1^{\prime}} [\mathbf{Q}=(π,π)]$", 
    title="B=$(Int(peierls)),  L=8", xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0,1)))

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("CDS_ππ")] , 
            dfs[!,Symbol("ΔCDS_ππ")]; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("CDS_ππ")]  ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10)
end
axislegend( position=(1, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_charge_A1p_ππ.png"), fig)


################
## A1' bilinear susceptibility
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{A_1^{\prime}} [\mathbf{Q}=(0,0)]$", 
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
 
    eb=errorbars!(dfs[!,:U], df_ϕs[!,Symbol("A1p_dQ_S_00")], 
            dfs[!,Symbol("ΔA1p_dQ_S_00")], df_ϕs[!,Symbol("ΔA1p_dQ_S_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:U], dfs[!,Symbol("A1p_dQ_S_00")]  ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10)

    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_S_00")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15)         
 
end

axislegend( position=(0, 1))
display(fig)
CairoMakie.save(joinpath(p, "chi_bilinear_A1p_00.png"), fig)

################
## A1' bilinear susceptibility as a function of T
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$χ_{\mathrm{bilinear}}^{A_1^{\prime}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(dfs[!,:T], df_ϕs[!,Symbol("A1p_dQ_S_00")]  , 
            dfs[!,Symbol("ΔA1p_dQ_S_00")], df_ϕs[!,Symbol("ΔA1p_dQ_S_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatter!(top, dfs[!,:T], dfs[!,Symbol("A1p_dQ_S_00")]  ;  marker = :xcross,
            color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10)

    CairoMakie.scatterlines!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("A1p_dQ_S_00")] ;  
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=15)            
              
end

axislegend( position=(1, 1))
display(fig)
CairoMakie.save(joinpath(p, "chi_bilinear_A1p_00_T.png"), fig)


################
## A1' bilinear correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$\mathrm{corr:\,}S_{\mathrm{bilinear}}^{A_1^{\prime}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
 
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("A1p_dQ_C_00")] , 
            dfs[!,Symbol("ΔA1p_dQ_C_00")], df_ϕs[!,Symbol("ΔA1p_dQ_C_00")]; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("A1p_dQ_C_00")]  ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_C_00")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15) 

end
axislegend( position=(0, 1))

################
## A1' bilinear order parameter 
################
top=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$\left| \Phi^{A_1^{\prime}} \right|$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_OP))  
 
    eb=errorbars!(df_OPs[!,:U], df_OPs[!,Symbol("A1p_OP")], 
    df_OPs[!,Symbol("ΔA1p_OP")], df_ϕ_OPs[!,Symbol("ΔA1p_OP")];
     linewidth=er_lw, whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:U], df_OPs[!,Symbol("A1p_OP")] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10) 
        
    CairoMakie.scatter!(top, df_ϕ_OPs[!,:U], df_ϕ_OPs[!,Symbol("A1p_OP")] ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15)
end
#axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "corr_bilinear_A1p_00.png"), fig)

################
## A1' bilinear correlation function as function of temperature
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], ylabel=L"$\mathrm{corr:\,}S_{\mathrm{bilinear}}^{A_1^{\prime}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30,xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=8:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("A1p_dQ_C_00")] , 
            dfs[!,Symbol("ΔA1p_dQ_C_00")], df_ϕs[!,Symbol("ΔA1p_dQ_C_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("A1p_dQ_C_00")]  ;  marker = :xcross,
            color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("A1p_dQ_C_00")] , 
            marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], 
             markersize=15)  
end
axislegend( position=(1, 1))

################
## A1' bilinear order parameter as function of temperature
################
top=Axis(fig[2, 1], xlabel=L"$T/t$", ylabel=L"$\left| \Phi^{A_1^{\prime}} \right|$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
for (idx, (Us, Ls)) in enumerate(zip(df_OP_LU[:,:U],df_OP_LU[:,:L]))
    num_T_points=length(df_OP_LU[:,:U]);
    df_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ_OP))  
 
    eb=errorbars!(df_OPs[!,:T], df_OPs[!,Symbol("A1p_OP")], 
        df_OPs[!,Symbol("ΔA1p_OP")], df_ϕ_OPs[!,Symbol("ΔA1p_OP")];
     linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:T], df_OPs[!,Symbol("A1p_OP")] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10)
        
    CairoMakie.scatter!(top, df_ϕ_OPs[!,:T], df_ϕ_OPs[!,Symbol("A1p_OP")] ;  
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=15) 
       
end
#axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "corr_bilinear_A1p_00_T.png"), fig)