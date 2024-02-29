


for U0 in [0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4]
# for U0 in [1.1, 2.0]


################
## Energy as a function of T
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$\mathrm{energy}$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),xticklabelsvisible=false,
    title=my_title
)

count=1;
for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    if Us ==U0
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  

    eb=errorbars!(top, dfs[!,:T], dfs[!,Symbol("energy")],
        dfs[!,Symbol("Δenergy")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("energy")] ;  
        marker = :xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "U=$(Us), L=$(Ls)"
    )
    count+=2;

    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("energy_TI")],
    dfs[!,Symbol("Δenergy_TI")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("energy_TI")] ;  
        marker = :xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "TI, U=$(Us), L=$(Ls)"
    )
    # count+=2;
    # CairoMakie.scatterlines!(top, rangeT, E_pot_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
    #     marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
    #     color = colorschemes[:tab20][count]
    # )
    end
end

axislegend( position=(1, 0))


################
## Specific heat as a function of T
################
bottom=Axis(fig[2, 1], xlabel=L"$T/t$", ylabel=L"$\mathrm{C_V}$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

    count=1;
for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    if Us ==U0
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  

    # eb=errorbars!(dfs[!,:T], dfs[!,Symbol("spec_heat")]  ,
    #     dfs[!,Symbol("Δspec_heat")] ; linewidth=er_lw, 
    #         whiskerwidth=10,color = colorschemes[:tab20][count])
    # CairoMakie.translate!(eb, 0, 0, -0.5)
    # CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("spec_heat")] ;  
    #     marker =:xcross, markersize=7,
    #     color = colorschemes[:tab20][count], label = "U=$(Us), L=$(Ls)")
    count+=2;
    eb=errorbars!(bottom, dfs[!,:T], dfs[!,Symbol("spec_heat_TI")]  ,
        dfs[!,Symbol("Δspec_heat_TI")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:T], dfs[!,Symbol("spec_heat_TI")] ;  
        marker =:xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "TI, U=$(Us), L=$(Ls)"
    )
    count+=2;
    CairoMakie.scatterlines!(bottom, rangeT, heat_capacity_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][count]
    )
    end
end

display(fig)
if save_bool_fig
    p_spec_heat=p * "spec_heat/"
    if !isdir(p_spec_heat)
        mkdir(p_spec_heat)
    end
    CairoMakie.save(joinpath(p_spec_heat, "spec_heat_U_" * to_string(U0) * ".png"), fig)
end





################
## h2
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$\mathrm{h2}$", xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),xticklabelsvisible=false,
    title=my_title2
)

count=1;
for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    if Us ==U0
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  

    eb=errorbars!(top, dfs[!,:T], dfs[!,Symbol("h2")] ./dfs[!,:T],
        dfs[!,Symbol("Δh2")] ./dfs[!,:T]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("h2")] ./dfs[!,:T] ;  
        marker = :xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "U=$(Us), L=$(Ls)"
    )
    count+=2;
    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("h2_TI")],
        dfs[!,Symbol("Δh2_TI")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("h2_TI")] ;  
        marker = :xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "TI, U=$(Us), L=$(Ls)"
    )

    count+=2;

    end
end

axislegend( position=(1, 1))


################
## h3
################
mid=Axis(fig[2, 1],  ylabel=L"$\mathrm{h3}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

count=1;
for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    if Us ==U0
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  

    eb=errorbars!(mid, dfs[!,:T], dfs[!,Symbol("h3")] ./dfs[!,:T],
        dfs[!,Symbol("Δh3")] ./dfs[!,:T] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, dfs[!,:T], dfs[!,Symbol("h3")] ./dfs[!,:T] ;  
        marker =:xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "U=$(Us), L=$(Ls)"
    )
    count+=2;
    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("h3_TI")], dfs[!,Symbol("Δh3_TI")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, dfs[!,:T], dfs[!,Symbol("h3_TI")] ;  
        marker =:xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "TI, U=$(Us), L=$(Ls)"
    )

    count+=2;

    end
end

################
## h4
################
bottom=Axis(fig[3, 1], xlabel=L"$T/t$", ylabel=L"$\mathrm{h4}$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing)
)

count=1;
for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    if Us ==U0
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(bottom, dfs[!,:T], dfs[!,Symbol("h4")] ,
        dfs[!,Symbol("Δh4")]  ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:T], dfs[!,Symbol("h4")] ;  
        marker =:xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "U=$(Us), L=$(Ls)"
    )
    count+=2;
    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("h4_TI")] ,
        dfs[!,Symbol("Δh4_TI")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:T], dfs[!,Symbol("h4_TI")] ;  
        marker =:xcross, markersize=7,
        color = colorschemes[:tab20][count], label = "TI, U=$(Us), L=$(Ls)"
    )

    count+=2;
    eb=errorbars!(df_ϕs[!,:T], df_ϕs[!,Symbol("h4ϕ")] ,
    df_ϕs[!,Symbol("Δh4ϕ")] ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][count]
    )
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_ϕs[!,:T], df_ϕs[!,Symbol("h4ϕ")] ;  
        marker =marker = '□', markersize=8,
        color = colorschemes[:tab20][count], label = "h4ϕ, U=$(Us), L=$(Ls)"
    )
    count+=2;
    CairoMakie.scatterlines!(bottom, rangeT, h4_analytic.(Us, rangeT .^(-1), μ0_plot, Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=1.5,
        color = colorschemes[:tab20][count]
    )
    end
end
axislegend( position=(1, 1))

display(fig)

if save_bool_fig
    p_hs=p * "h2_h3_h4/"
    if !isdir(p_hs)
        mkdir(p_hs)
    end
    CairoMakie.save(joinpath(p_hs, "h2_h3_h4_U_" * to_string(U0) * ".png"), fig)
end

end