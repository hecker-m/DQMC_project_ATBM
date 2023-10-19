df_Binder=sort(df_Binder, [:U, :T, :B, :L]);
df_Binder_LU=sort(unique(df_Binder[:,[:U,:L]]),:U)

################
## Magnetic Binder cumulant as a function of T
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$\mathrm{magn.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_Binder_LU[:,:U],df_Binder_LU[:,:L]))
    num_T_points=length(df_Binder_LU[:,:U]);
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_Binder))  

    eb=errorbars!(dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_magnetic")],
            dfs_Binder[!,Symbol("ΔBinder_magnetic")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_magnetic")] ;  
        marker = '□', markersize=15,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)")
         
end

axislegend( position=(1, 1))
display(fig)
CairoMakie.save(joinpath(p, "Binder_Magn_T.png"), fig)

################
## (Selected) Magnetic Binder cumulant as a function of T
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$\mathrm{magn.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0, 0.7)),
    title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_Binder_LU[:,:U],df_Binder_LU[:,:L]))
    num_T_points=length(df_Binder_LU[:,:U]);
    if Us >0.9
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_Binder))  

    eb=errorbars!(dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_magnetic")],
            dfs_Binder[!,Symbol("ΔBinder_magnetic")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_magnetic")] ;  
        marker = '□', markersize=15,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)")
    end
end

axislegend( position=(1, 1))
display(fig)
CairoMakie.save(joinpath(p, "Binder_Magn_T_select.png"), fig)



################
## Nematic Binder cumulant as a function of T
################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$T/t$", ylabel=L"$\mathrm{nem.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-3,3)),
    title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_Binder_LU[:,:U],df_Binder_LU[:,:L]))
    if Us >0.9
    num_T_points=length(df_Binder_LU[:,:U]);
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_Binder))  

    eb=errorbars!(dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_nematic")],
            dfs_Binder[!,Symbol("ΔBinder_nematic")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_nematic")] ;  
        marker = '□', markersize=15,
        color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)")
    end
end

axislegend( position=(1, 1))
display(fig)
CairoMakie.save(joinpath(p, "Binder_Nematic_T.png"), fig)