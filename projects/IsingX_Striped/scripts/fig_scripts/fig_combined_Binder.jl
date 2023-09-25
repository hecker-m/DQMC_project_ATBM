include("../../../../src/MonteCarlo.jl_modified/src/my_files/helpful_fcns.jl")

df_Binder=sort(df_Binder, [:U, :T, :B, :L]);
df_Binder_LU=sort(unique(df_Binder[:,[:U,:L]]),:U)

for U0 in [1.0, 1.1, 1.3, 1.6, 2.0]


################
## Magnetic Binder cumulant as a function of T
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$\mathrm{magn.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0.3, nothing)),xticklabelsvisible=false,
    title=L"B=%$(Int(peierls)) \;:\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [xx/10 for xx in 1:2:9], color = :gray, linewidth=0.2), 0, 0, -0.8)
# CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)
count=1;
for (idx, (Us, Ls)) in enumerate(zip(df_Binder_LU[:,:U],df_Binder_LU[:,:L]))
    if Us ==U0
    num_T_points=length(df_Binder_LU[:,:U]);
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_Binder))  

    eb=errorbars!(dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_magnetic")],
            dfs_Binder[!,Symbol("ΔBinder_magnetic")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_magnetic")] ;  
        marker = '□', markersize=15,
        color = colorschemes[:tab20][count], label = "U=$(Us), L=$(Ls)")
    count+=2;
    end
end

axislegend( position=(1, 1))

################
## Nematic Binder cumulant as a function of T
################
top=Axis(fig[2, 1], xlabel=L"$T/t$", ylabel=L"$\mathrm{nem.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-1.5,1)))

    CairoMakie.translate!(vlines!(top, [xx/10 for xx in 1:2:9], color = :gray, linewidth=0.2), 0, 0, -0.8)
    #CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

    count=1;
for (idx, (Us, Ls)) in enumerate(zip(df_Binder_LU[:,:U],df_Binder_LU[:,:L]))
    if Us ==U0
    num_T_points=length(df_Binder_LU[:,:U]);
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_Binder))  

    eb=errorbars!(dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_nematic")],
            dfs_Binder[!,Symbol("ΔBinder_nematic")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:T], dfs_Binder[!,Symbol("Binder_nematic")] ;  
        marker = '□', markersize=15,
        color = colorschemes[:tab20][count], label = "U=$(Us), L=$(Ls)")
    count+=2;
    end
end

display(fig)
CairoMakie.save(joinpath(p, "Binder_Combined_U_" * to_string(U0) * ".png"), fig)


end