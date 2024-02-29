   
path_Analysis_Fcns="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/Analysis_Fcns/"
include(path_Analysis_Fcns * "analytical_fcns.jl")


df_Binder=sort(df_Binder, [:U, :T, :B, :L]);
df_Binder_LU=sort(unique(df_Binder[:,[:U,:L]]),:U)
df_Binder_LT=sort(unique(df_Binder[:,[:T,:L]]),:T, rev=true)

β_list=1.0 ./ unique(df[!,:T])   
for β0 in β_list

    ################
    ## h4
    ################
    fig = Figure(resolution = (800, 800))

    top=Axis(fig[1, 1], title=mytitle, ylabel=L"$\mathrm{h}_4/U^2$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, (2, 15)))
    mid=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$Δ\mathrm{h}_4 /U^2$", 
        xlabelsize=30, ylabelsize=30, 
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.16, 0.16)))  
        
        count=1;
    for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
        if isapprox(Ts, 1/β0)
        dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
        df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
    
        minU=findmin(df[!,:U])[1];
        maxU=findmax(df[!,:U])[1];
        rangeU=range(minU, maxU, 30)
        CairoMakie.lines!(mid, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

        eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("h4")] ./(dfs[!,:U].^2),
            dfs[!,Symbol("Δh4")]  ./(dfs[!,:U].^2); linewidth=er_lw, 
                whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(eb, 0, 0, -0.5)
        CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("h4")]  ./(dfs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")

        eb=errorbars!(mid, dfs[!,:U], dfs[!,Symbol("h4")]  ./(dfs[!,:U].^2) .-
            h4_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2),
            dfs[!,Symbol("Δh4")]  ./(dfs[!,:U].^2); linewidth=er_lw, 
                whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(eb, 0, 0, -0.5)
        CairoMakie.scatterlines!(mid, dfs[!,:U], dfs[!,Symbol("h4")]  ./(dfs[!,:U].^2) .-
            h4_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")   

        count+=2;
        eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("h4_TI")] ./(dfs[!,:U].^2),
        dfs[!,Symbol("Δh4_TI")] ./(dfs[!,:U].^2); linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(eb, 0, 0, -0.5)
        CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("h4_TI")] ./(dfs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$(TI):\;β=%$(β0), L=%$(Ls) $")

        eb=errorbars!(mid, dfs[!,:U], dfs[!,Symbol("h4_TI")] ./(dfs[!,:U].^2) .-
            h4_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2),
            dfs[!,Symbol("Δh4_TI")] ./(dfs[!,:U].^2); linewidth=er_lw, 
                whiskerwidth=10,color = colorschemes[:tab20][count])
            CairoMakie.translate!(eb, 0, 0, -0.5)            
        CairoMakie.scatterlines!(mid, dfs[!,:U], dfs[!,Symbol("h4_TI")] ./(dfs[!,:U].^2) .-
            h4_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$(TI):\;β=%$(β0), L=%$(Ls) $") 
            
        count+=2;
        eb=errorbars!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("h4ϕ")] ./(df_ϕs[!,:U].^2),
        df_ϕs[!,Symbol("Δh4ϕ")] ./(df_ϕs[!,:U].^2); linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(eb, 0, 0, -0.5)
        CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("h4ϕ")] ./(df_ϕs[!,:U].^2);  
        marker =:xcross, markersize=8,
        color = colorschemes[:tab20][count], label = L"$(ϕ):\;β=%$(β0), L=%$(Ls) $")

        eb=errorbars!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("h4ϕ")] ./(df_ϕs[!,:U].^2) .-
            h4_analytic.(df_ϕs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(df_ϕs[!,:U].^2),
            df_ϕs[!,Symbol("Δh4ϕ")] ./(df_ϕs[!,:U].^2); linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(eb, 0, 0, -0.5)
        CairoMakie.scatterlines!(mid, df_ϕs[!,:U],  df_ϕs[!,Symbol("h4ϕ")] ./(df_ϕs[!,:U].^2) .-
            h4_analytic.(df_ϕs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(df_ϕs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$(ϕ):\;β=%$(β0), L=%$(Ls) $") 

        count+=2;
        CairoMakie.scatterlines!(top, rangeU, h4_analytic.(rangeU, β0, 0.5*(μ0-(-μ0)), Ls) ./(rangeU.^2);  
            marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
            color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
       
        end
    end
    axislegend(top, position=(0, 1))
    axislegend(mid, position=(0.2, 0))

    display(fig)
    if save_bool_fig
        CairoMakie.save(joinpath(p, "h4_t0_" * to_string(δτ) * ".png"), fig)
    end


    ################
    ## h4 onsite
    ################
    fig = Figure(resolution = (800, 800))

    bottom=Axis(fig[1, 1], title=mytitle, ylabel=L"$\mathrm{h}_4^{\mathrm{OnSite}}/U^2$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing,nothing))   
    truebottom=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$Δ\mathrm{h}_4^{\mathrm{OnSite}}/U^2$", 
        xlabelsize=30, ylabelsize=30,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing,(-0.01, 0.01)))

        count=1;
    for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
        if isapprox(Ts, 1/β0)
        dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
        df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  
    
        minU=findmin(df[!,:U])[1];
        maxU=findmax(df[!,:U])[1];
        rangeU=range(minU, maxU, 30)
        CairoMakie.lines!(truebottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

        ebOS=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("h4_onsite")]  ./(dfs[!,:U].^2),
            dfs[!,Symbol("Δh4_onsite")]  ./(dfs[!,:U].^2); linewidth=er_lw, 
                whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(ebOS, 0, 0, -0.5)
        CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("h4_onsite")]  ./(dfs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")

        ebOS=errorbars!(truebottom, dfs[!,:U], dfs[!,Symbol("h4_onsite")]  ./(dfs[!,:U].^2).-
            h4_onsite_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2),
            dfs[!,Symbol("Δh4_onsite")]  ./(dfs[!,:U].^2); linewidth=er_lw, 
                whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(ebOS, 0, 0, -0.5)
        CairoMakie.scatterlines!(truebottom, dfs[!,:U], dfs[!,Symbol("h4_onsite")]  ./(dfs[!,:U].^2) .-
            h4_onsite_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")   

        count+=2;

        ebOSTI=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("h4_onsite_TI")] ./(dfs[!,:U].^2),
                dfs[!,Symbol("Δh4_onsite_TI")] ./(dfs[!,:U].^2); linewidth=er_lw, 
                    whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(ebOSTI, 0, 0, -0.5)
        CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("h4_onsite_TI")] ./(dfs[!,:U].^2);  
        marker =:xcross, markersize=7,
        color = colorschemes[:tab20][count], label = L"$(TI):\;β=%$(β0), L=%$(Ls) $")

        ebOSTI=errorbars!(truebottom, dfs[!,:U], dfs[!,Symbol("h4_onsite_TI")] ./(dfs[!,:U].^2) .-
            h4_onsite_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2),
            dfs[!,Symbol("Δh4_onsite_TI")] ./(dfs[!,:U].^2); linewidth=er_lw, 
                whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(ebOSTI, 0, 0, -0.5)
        CairoMakie.scatterlines!(truebottom, dfs[!,:U], dfs[!,Symbol("h4_onsite_TI")] ./(dfs[!,:U].^2) .-
            h4_onsite_analytic.(dfs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(dfs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$(TI):\;β=%$(β0), L=%$(Ls) $") 
            
        count+=2;

        ebOSϕ=errorbars!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("h4ϕ_OnSite")] ./(df_ϕs[!,:U].^2),
                df_ϕs[!,Symbol("Δh4ϕ_OnSite")] ./(df_ϕs[!,:U].^2); linewidth=er_lw, 
                    whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(ebOSϕ, 0, 0, -0.5)
        CairoMakie.scatterlines!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("h4ϕ_OnSite")] ./(df_ϕs[!,:U].^2);  
        marker =marker = '□', markersize=8,
        color = colorschemes[:tab20][count], label = L"$(ϕ):\;β=%$(β0), L=%$(Ls) $")

        ebOSϕ=errorbars!(truebottom, df_ϕs[!,:U], df_ϕs[!,Symbol("h4ϕ_OnSite")] ./(df_ϕs[!,:U].^2) .-
            h4_onsite_analytic.(df_ϕs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(df_ϕs[!,:U].^2),
            df_ϕs[!,Symbol("Δh4ϕ_OnSite")] ./(df_ϕs[!,:U].^2); linewidth=er_lw, 
                whiskerwidth=10,color = colorschemes[:tab20][count])
        CairoMakie.translate!(ebOSϕ, 0, 0, -0.5)
        CairoMakie.scatterlines!(truebottom, dfs[!,:U],  df_ϕs[!,Symbol("h4ϕ_OnSite")] ./(df_ϕs[!,:U].^2) .-
            h4_onsite_analytic.(df_ϕs[!,:U], β0, 0.5*(μ0-(-μ0)), Ls) ./(df_ϕs[!,:U].^2);  
            marker =:xcross, markersize=7,
            color = colorschemes[:tab20][count], label = L"$(ϕ):\;β=%$(β0), L=%$(Ls) $") 

        count+=2;
        CairoMakie.scatterlines!(bottom, rangeU, h4_onsite_analytic.(rangeU, β0, 0.5*(μ0-(-μ0)), Ls) ./(rangeU.^2);  
            marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
            color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")  
      
        end
    end
    axislegend(bottom, position=(1, 0))
    axislegend(truebottom, position=(1, 0))

    display(fig)
    if save_bool_fig
        CairoMakie.save(joinpath(p, "h4_onsite_t0_" * to_string(δτ) * ".png"), fig)
    end
end

################
## spin correlation function, S_A1 ∼ Q=(π,0) + Q=(0,π)
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$S_{\mathrm{spin}}^{A_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0.04, 0.14)))
mid=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$ΔS_{\mathrm{spin}}^{A_1}$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.01, 0.01)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(mid, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    eb=errorbars!(top, dfs[!,:U], (dfs[!,Symbol("SDC_A1_Mx_x")] ), 
        (dfs[!,Symbol("ΔSDC_A1_Mx_x")] )  ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], (dfs[!,Symbol("SDC_A1_Mx_x")])  ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    eb=errorbars!(mid, dfs[!,:U], (dfs[!,Symbol("SDC_A1_Mx_x")] )  .-
        Sspin_A1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        (dfs[!,Symbol("ΔSDC_A1_Mx_x")] ) ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, dfs[!,:U], (dfs[!,Symbol("SDC_A1_Mx_x")])   .-
        Sspin_A1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    eb=errorbars!(top, dfs[!,:U], df_ϕs[!,Symbol("SDC_A1_00")],  
        df_ϕs[!,Symbol("ΔSDC_A1_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("SDC_A1_00")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    eb=errorbars!(mid, dfs[!,:U], df_ϕs[!,Symbol("SDC_A1_00")] .-
            Sspin_A1_XX_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls),  
            df_ϕs[!,Symbol("ΔSDC_A1_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("SDC_A1_00")] .-
            Sspin_A1_XX_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, Sspin_A1_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $") 
end
end
axislegend(top, position=(1, 0))
axislegend(mid, position=(0, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Sspin_A1_t0_" * to_string(δτ) * ".png"), fig)
end

################
## spin correlation (2) function, S(2)_A1 
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$S_{\mathrm{spin}}^{(2),A_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0,0.035)))
mid=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$ΔS_{\mathrm{spin}}^{(2),A_1}$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.01,0.01)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df_ϕs[!,:U])[1];
    maxU=findmax(df_ϕs[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(mid, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 


    eb=errorbars!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("Sspin_2_A1_xx")],  
        df_ϕs[!,Symbol("ΔSspin_2_A1_xx")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("Sspin_2_A1_xx")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    eb=errorbars!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("Sspin_2_A1_xx")] .-
            Sspin2_A1_XX_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls),  
            df_ϕs[!,Symbol("ΔSspin_2_A1_xx")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("Sspin_2_A1_xx")] .-
            Sspin2_A1_XX_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, Sspin2_A1_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $") 
end
end
axislegend(top, position=(1, 0))
axislegend(mid, position=(0, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Sspin2_A1_t0_" * to_string(δτ) * ".png"), fig)
end

################
## magnetic Binder cumulant
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$\mathrm{magn.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0, 1)))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$ΔU_L$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.2, 0.2)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_Binder_LT[:,:T],df_Binder_LT[:,:L]))
    num_U_points=length(df_Binder_LT[:,:T]);
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_Binder)) 

    minU=findmin(df_Binder[!,:U])[1];
    maxU=findmax(df_Binder[!,:U])[1];
    minU=0.01;
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_magnetic")] , 
        dfs_Binder[!,Symbol("ΔBinder_magnetic")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_magnetic")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_magnetic")] .-
        UL_magnetic_analytic.(dfs_Binder[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs_Binder[!,Symbol("ΔBinder_magnetic")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_magnetic")] .-
        UL_magnetic_analytic.(dfs_Binder[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, UL_magnetic_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(1, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Binder_magnetic_t0_" * to_string(δτ) * ".png"), fig)
end

################
## spin susceptibility, χ_A1 ∼ Q=(π,0) + Q=(0,π)
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$χ_{\mathrm{spin}}^{A_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
mid=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δχ_{\mathrm{spin}}^{A_1}$", 
    xlabelsize=30, ylabelsize=30,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(mid, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    eb=errorbars!(top, dfs[!,:U], (dfs[!,Symbol("SDS_A1_Mx_x")] ) , 
        (dfs[!,Symbol("ΔSDS_A1_Mx_x")] )  ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], (dfs[!,Symbol("SDS_A1_Mx_x")])   ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    eb=errorbars!(mid, dfs[!,:U], (dfs[!,Symbol("SDS_A1_Mx_x")] ) .-
        χspin_A1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        (dfs[!,Symbol("ΔSDS_A1_Mx_x")] )  ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, dfs[!,:U], (dfs[!,Symbol("SDS_A1_Mx_x")] )   .-
        χspin_A1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    eb=errorbars!(top, dfs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")],  
        df_ϕs[!,Symbol("ΔSDS_A1_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    eb=errorbars!(mid, dfs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")] .-
        χspin_A1_XX_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls),  
            df_ϕs[!,Symbol("ΔSDS_A1_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")] .-
            χspin_A1_XX_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χspin_A1_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $") 
end
end
axislegend(top, position=(1, 0))
axislegend(mid, position=(0, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_spin_A1_t0_" * to_string(δτ) * ".png"), fig)
end

################
## potential energy
################
fig = Figure(resolution = (800, 800))
mid=Axis(fig[1, 1], title=mytitle,   ylabel=L"$-E_{\mathrm{pot}}/U$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$-ΔE_{\mathrm{pot}}/U$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.03, 0.03)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(mid, dfs[!,:U], -dfs[!,Symbol("pot_energy")] ./dfs[!,:U] /(Ls^2), 
        dfs[!,Symbol("Δpot_energy")]./dfs[!,:U] /(Ls^2);
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, dfs[!,:U], -dfs[!,Symbol("pot_energy")] ./dfs[!,:U] /(Ls^2);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)


    ebbot=errorbars!(bottom, dfs[!,:U], -dfs[!,Symbol("pot_energy")] ./dfs[!,:U] /(Ls^2) .+
        E_pot_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ./dfs[!,:U], 
        dfs[!,Symbol("Δpot_energy")]./dfs[!,:U] /(Ls^2);
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], -dfs[!,Symbol("pot_energy")] ./dfs[!,:U] /(Ls^2) .+
        E_pot_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ./dfs[!,:U];  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    ebbot=errorbars!(mid, df_ϕs[!,:U], -df_ϕs[!,Symbol("E_potϕ")] ./df_ϕs[!,:U], 
            df_ϕs[!,Symbol("ΔE_potϕ")]./df_ϕs[!,:U];
             linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(mid, df_ϕs[!,:U], -df_ϕs[!,Symbol("E_potϕ")] ./df_ϕs[!,:U] ;  
            marker = :xcross, color = colorschemes[:tab20][count], 
            label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, df_ϕs[!,:U], -df_ϕs[!,Symbol("E_potϕ")] ./df_ϕs[!,:U] .+
        E_pot_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ./df_ϕs[!,:U], 
        df_ϕs[!,Symbol("ΔE_potϕ")]./df_ϕs[!,:U] ;
            linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_ϕs[!,:U], -df_ϕs[!,Symbol("E_potϕ")] ./df_ϕs[!,:U] .+
        E_pot_analytic.(df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ./df_ϕs[!,:U];  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(mid, rangeU, -E_pot_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ./(rangeU);  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $") 

end
end
axislegend(mid, position=(1, 0))
axislegend(bottom, position=(0.2, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "E_pot_t0_" * to_string(δτ) * ".png"), fig)
end


################
## B1 nematic correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$\mathrm{corr:\,}S_{\mathrm{bilinear}}^{B_1} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0, 0.018)))
bottom=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$ΔS_{\mathrm{bilinear}}^{B_1}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.005, 0.005)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("NemC_X")] , dfs[!,Symbol("ΔNemC_X")]; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("NemC_X")]  ;  marker = :xcross,
        color = colorschemes[:tab20][count], markersize=10,
        label=L"$():\;β=%$(1/Ts), L=%$(Ls) $")

    eb=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("NemC_X")] .-
    Sbil_B1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), dfs[!,Symbol("ΔNemC_X")]; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("NemC_X")] .-
        Sbil_B1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  marker = :xcross,
        color = colorschemes[:tab20][count], markersize=10,
        label=L"$():\;β=%$(1/Ts), L=%$(Ls) $")

    count+=2;
    ebϕ=errorbars!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("NemC")] ,df_ϕs[!,Symbol("ΔNemC")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("NemC")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    ebϕ=errorbars!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("NemC")] .- 
    Sbil_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) , df_ϕs[!,Symbol("ΔNemC")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("NemC")] .- 
        Sbil_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, Sbil_B1_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")

end
axislegend(top, position=(1, 0))
axislegend(bottom, position=(0.3, 1))

end
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Sbil_B1_t0_" * to_string(δτ) * ".png"), fig)
end

################
## Bilinear (2) B1 correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$S_{\mathrm{bilinear}}^{(2),B_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.005,0.005)))
bottom=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$ΔS_{\mathrm{bilinear}}^{(2),B_1}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.002,0.002)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df_ϕs[!,:U])[1];
    maxU=findmax(df_ϕs[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebϕ=errorbars!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_B1_xx")] ,df_ϕs[!,Symbol("ΔSbil_2_B1_xx")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_B1_xx")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    ebϕ=errorbars!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_B1_xx")] .- 
        Sbil2_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) , df_ϕs[!,Symbol("ΔSbil_2_B1_xx")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_B1_xx")] .- 
        Sbil2_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, Sbil2_B1_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")

end
axislegend(top, position=(1, 0))
axislegend(bottom, position=(0.3, 1))

end
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Sbil2_B1_t0_" * to_string(δτ) * ".png"), fig)
end

################
## B1 bilinear Binder cumulant
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$B_1 \mathrm{-bil.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-3,1)))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$ΔU_L$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-1,1)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_Binder_LT[:,:T],df_Binder_LT[:,:L]))
    num_U_points=length(df_Binder_LT[:,:T]);
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_Binder)) 

    minU=findmin(df_Binder[!,:U])[1];
    maxU=findmax(df_Binder[!,:U])[1];
    minU=0.01;
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_nematic")] , 
        dfs_Binder[!,Symbol("ΔBinder_nematic")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_nematic")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_nematic")] .-
        UL_bil_B1_analytic.(dfs_Binder[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs_Binder[!,Symbol("ΔBinder_nematic")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_nematic")] .-
        UL_bil_B1_analytic.(dfs_Binder[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, UL_bil_B1_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(1, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Binder_bil_B1_t0_" * to_string(δτ) * ".png"), fig)
end


################
## B1 nematic susceptibility
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$χ_{\mathrm{bilinear}}^{B_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0, 0.018)))
bottom=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$Δχ_{\mathrm{bilinear}}^{B_1}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.003, 0.003)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("NemS_X")] , dfs[!,Symbol("ΔNemS_X")]; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("NemS_X")]  ;  marker = :xcross,
        color = colorschemes[:tab20][count], markersize=10,
        label=L"$():\;β=%$(1/Ts), L=%$(Ls) $")

    eb=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("NemS_X")] .-
        χbil_B1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), dfs[!,Symbol("ΔNemS_X")]; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("NemS_X")] .-
        χbil_B1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  marker = :xcross,
        color = colorschemes[:tab20][count], markersize=10,
        label=L"$():\;β=%$(1/Ts), L=%$(Ls) $")

    count+=2;
    ebϕ=errorbars!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("NemS")] ,df_ϕs[!,Symbol("ΔNemS")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("NemS")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    ebϕ=errorbars!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("NemS")] .- 
        χbil_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) , df_ϕs[!,Symbol("ΔNemS")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("NemS")] .- 
        χbil_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χbil_B1_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")

end
axislegend(top, position=(1, 0))
axislegend(bottom, position=(0.3, 1))

end
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_bil_B1_t0_" * to_string(δτ) * ".png"), fig)
end


################
## A1` bilinear correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$\mathrm{corr:\,}S_{\mathrm{bilinear}}^{A_1^{′}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (0, 0.018)))
bottom=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$ΔS_{\mathrm{bilinear}}^{A_1^{′}}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.005, 0.005)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    eb=errorbars!(top, dfs[!,:U], dfs[!,Symbol("A1p_dQ_C_X")] , dfs[!,Symbol("ΔA1p_dQ_C_X")]; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("A1p_dQ_C_X")]  ;  marker = :xcross,
        color = colorschemes[:tab20][count], markersize=10,
        label=L"$():\;β=%$(1/Ts), L=%$(Ls) $")

    eb=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("A1p_dQ_C_X")] .-
    Sbil_B1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), dfs[!,Symbol("ΔA1p_dQ_C_X")]; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("A1p_dQ_C_X")] .-
        Sbil_B1_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  marker = :xcross,
        color = colorschemes[:tab20][count], markersize=10,
        label=L"$():\;β=%$(1/Ts), L=%$(Ls) $")

    count+=2;
    ebϕ=errorbars!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_C")] ,df_ϕs[!,Symbol("ΔA1p_dQ_C")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_C")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    ebϕ=errorbars!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_C")] .- 
    Sbil_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) , df_ϕs[!,Symbol("ΔA1p_dQ_C")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("A1p_dQ_C")] .- 
        Sbil_B1_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, Sbil_A1P_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")

end
axislegend(top, position=(1, 0))
axislegend(bottom, position=(0.3, 1))

end
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Sbil_A1P_t0_" * to_string(δτ) * ".png"), fig)
end

################
## Bilinear (2) A1` correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,  ylabel=L"$S_{\mathrm{bilinear}}^{(2),A_1^{\prime}}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.005,0.005)))
bottom=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$ΔS_{\mathrm{bilinear}}^{(2),A_1^{\prime}}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.002,0.002)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df_ϕs[!,:U])[1];
    maxU=findmax(df_ϕs[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebϕ=errorbars!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] ,df_ϕs[!,Symbol("ΔSbil_2_A1p_xx")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] ;  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    ebϕ=errorbars!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] .- 
        Sbil2_A1P_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) , df_ϕs[!,Symbol("ΔSbil_2_A1p_xx")] ; 
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebϕ, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, df_ϕs[!,:U], df_ϕs[!,Symbol("Sbil_2_A1p_xx")] .- 
        Sbil2_A1P_XX_analytic.( df_ϕs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
            marker = :xcross, color = colorschemes[:tab20][count],  markersize=10,
            label=L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $") 

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, Sbil2_A1P_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")

end
axislegend(top, position=(1, 0))
axislegend(bottom, position=(0.3, 1))

end
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Sbil2_A1P_t0_" * to_string(δτ) * ".png"), fig)
end

################
## A1` bilinear Binder cumulant
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$A^{\prime}_1 \mathrm{-bil.\,Binder}\;U_L$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-3,1)))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$ΔU_L$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-1,1)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_Binder_LT[:,:T],df_Binder_LT[:,:L]))
    num_U_points=length(df_Binder_LT[:,:T]);
    dfs_Binder=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_Binder)) 

    minU=findmin(df_Binder[!,:U])[1];
    maxU=findmax(df_Binder[!,:U])[1];
    minU=0.01;
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_A1p")] , 
        dfs_Binder[!,Symbol("ΔBinder_A1p")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_A1p")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_A1p")] .-
        UL_bil_A1P_analytic.(dfs_Binder[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs_Binder[!,Symbol("ΔBinder_A1p")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs_Binder[!,:U], dfs_Binder[!,Symbol("Binder_A1p")] .-
        UL_bil_A1P_analytic.(dfs_Binder[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(ϕ):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, UL_bil_A1P_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(1, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Binder_bil_A1P_t0_" * to_string(δτ) * ".png"), fig)
end



################
## compressibility
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$χ_{\mathrm{charge}} $", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δχ_{\mathrm{charge}} $", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    maxU=findmax(df[!,:U])[1];
    minU=0.001;
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("compress")] , 
        dfs[!,Symbol("Δcompress")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("compress")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("compress")] .-
        χ_charge_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("Δcompress")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("compress")] .-
        χ_charge_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χ_charge_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(1, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Chi_charge_t0_" * to_string(δτ) * ".png"), fig)
end

################
## χ_pair
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$χ_{\mathrm{pair}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δχ_{\mathrm{pair}} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.04,0.04)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    minU=0.001
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("PDS_spm_00")] , 
        dfs[!,Symbol("ΔPDS_spm_00")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("PDS_spm_00")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(zz):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("PDS_spm_00")] .-
        χ_pair_ZZ_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("ΔPDS_spm_00")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("PDS_spm_00")] .-
        χ_pair_ZZ_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(zz):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χ_pair_ZZ_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(zz,ana):\;β=%$(1/Ts), L=%$(Ls) $")

    count+=2;
    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("PDS_s_00")] , 
        dfs[!,Symbol("ΔPDS_s_00")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("PDS_s_00")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(00):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("PDS_s_00")] .-
        χ_pair_00_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("ΔPDS_s_00")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("PDS_s_00")] .-
        χ_pair_00_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$(00):\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χ_pair_00_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(00,ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(1, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Chi_pair_t0_" * to_string(δτ) * ".png"), fig)
end



################
## specific heat
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$\mathrm{C_V}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δ\mathrm{C_V}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    minU=0.01;
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("spec_heat_TI")] , 
        dfs[!,Symbol("Δspec_heat_TI")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("spec_heat_TI")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("spec_heat_TI")] .-
        heat_capacity_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("Δspec_heat_TI")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("spec_heat_TI")] .-
        heat_capacity_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, heat_capacity_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(0, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "heat_capacity_t0_" * to_string(δτ) * ".png"), fig)
end




################
## phase stiffness
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$\rho_{s}\;\cdot  β\,\frac{\pi}{2} $", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δ\rho_{s}\;\cdot  β\,\frac{\pi}{2} $", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, (-0.0004, 0.0004)))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    minU=0.01;
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("ρ_s")] *π/(2Ts), 
        dfs[!,Symbol("Δρ_s")] *π/(2Ts);
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("ρ_s")] *π/(2Ts);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("ρ_s")] .*π/(2Ts) .-
        heat_capacity_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("Δρ_s")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("ρ_s")] .* π/(2Ts) .-
        heat_capacity_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)


end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(0, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "phase_stiffness_t0_" * to_string(δτ) * ".png"), fig)
end


################
## charge-B1 susceptibility
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$χ_{\mathrm{charge}}^{B_1} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δχ_{\mathrm{charge}}^{B_1} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    minU=0.01;
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("B1_CDS_00")] .*Ts, 
        dfs[!,Symbol("ΔB1_CDS_00")] .*Ts;
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("B1_CDS_00")] .*Ts;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("B1_CDS_00")] .*Ts .-
        χ_charge_B1_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("ΔB1_CDS_00")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("B1_CDS_00")] .*Ts .-
        χ_charge_B1_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χ_charge_B1_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(0, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_charge_B1_t0_" * to_string(δτ) * ".png"), fig)
end

################
## charge-B1 correlation function
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$S_{\mathrm{charge}}^{B_1} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$ΔS_{\mathrm{charge}}^{B_1} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    minU=0.01;
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("B1_CDC_00")] , 
        dfs[!,Symbol("ΔB1_CDC_00")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("B1_CDC_00")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("B1_CDC_00")] .-
        S_charge_B1_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("ΔB1_CDC_00")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("B1_CDC_00")] .-
        S_charge_B1_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, S_charge_B1_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(0, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "corr_charge_B1_t0_" * to_string(δτ) * ".png"), fig)
end


################
## charge-A1P susceptibility
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$χ_{\mathrm{charge}}^{A_1^{′}}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δχ_{\mathrm{charge}}^{A_1^{′}}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    minU=0.01;
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("CDS_ππ")] , 
        dfs[!,Symbol("ΔCDS_ππ")];
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("CDS_ππ")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("CDS_ππ")] .-
        χ_charge_A1P_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) *Ls^2, 
        dfs[!,Symbol("ΔCDS_ππ")];
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("CDS_ππ")] .-
        χ_charge_A1P_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls) *Ls^2;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χ_charge_A1P_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) *Ls^2 ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")
end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(1, 0))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_charge_A1P_t0_" * to_string(δτ) * ".png"), fig)
end




################
## charge-XX susceptibility
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], title=mytitle,   ylabel=L"$χ_{\mathrm{charge}}^{XX} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$U/t$", ylabel=L"$Δχ_{\mathrm{charge}}^{XX} [\mathbf{Q}=(0,0)]$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1, ]
    count=1;
for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    minU=findmin(df[!,:U])[1];
    minU=0.01;
    maxU=findmax(df[!,:U])[1];
    rangeU=range(minU, maxU, 30)
    CairoMakie.lines!(bottom, rangeU, zeros(length(rangeU)) ; linewidth=1.5, color = :black) 

    ebbot=errorbars!(top, dfs[!,:U], dfs[!,Symbol("CDSxx_00")] *Ts/Ls^2, 
        dfs[!,Symbol("ΔCDSxx_00")] *Ts/Ls^2,
         linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("CDSxx_00")] *Ts/Ls^2;  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    ebbot=errorbars!(bottom, dfs[!,:U], dfs[!,Symbol("CDSxx_00")] *Ts/Ls^2 .-
        χ_charge_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls), 
        dfs[!,Symbol("ΔCDSxx_00")] *Ts/Ls^2;
        linewidth=er_lw,   whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(ebbot, 0, 0, -0.5)
    CairoMakie.scatterlines!(bottom, dfs[!,:U], dfs[!,Symbol("CDSxx_00")] *Ts/Ls^2 .-
        χ_charge_XX_analytic.(dfs[!,:U], 1/Ts, 0.5*(μ0-(-μ0)), Ls);  
        marker = :xcross, color = colorschemes[:tab20][count], 
        label = L"$():\;β=%$(1/Ts), L=%$(Ls) $", markersize=10)

    count+=2;
    CairoMakie.scatterlines!(top, rangeU, χ_charge_XX_analytic.(rangeU, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(1/Ts), L=%$(Ls) $")


end
end
axislegend(top, position=(1, 1))
axislegend(bottom, position=(0, 1))
display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "chi_charge_XX_t0_" * to_string(δτ) * ".png"), fig)
end