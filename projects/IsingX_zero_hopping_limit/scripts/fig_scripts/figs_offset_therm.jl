   
path_Analysis_Fcns="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/Analysis_Fcns/"
include(path_Analysis_Fcns * "analytical_fcns.jl")


gdf=groupby(df, :δτ)
gdf_ϕ=groupby(df_ϕ, :δτ)

ms=3;lw=0.8;
################
## potential energy
################
for dummy in [1,]

    fig = Figure(resolution = (800, 800))


    mid=Axis(fig[1, 1], title=mytitle,   ylabel=L"$E_{\mathrm{pot}}$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    bottom=Axis(fig[2, 1],  xlabel=L"$sweep*10$", ylabel=L"$E^{ϕ}_{\mathrm{pot}}$", 
        xlabelsize=30, ylabelsize=30, 
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    
    count=1;
    for n in 1:length(gdf)
        Ts=gdf[n][1, :T]
        Us=gdf[n][1, :U]
        Ls=gdf[n][1, :L]
        δτ=gdf[n][1, :δτ]

        CairoMakie.scatterlines!(mid, gdf[n][!,:_th_meas], gdf[n][!,Symbol("pot_energy")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)

        CairoMakie.scatterlines!(bottom, gdf_ϕ[n][!,:_th_meas], gdf_ϕ[n][!,Symbol("pot_energy_ϕ")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)



        if count==1
            range_sweeps=1:200;
            range_sweeps=range(1, 200, 30)
            CairoMakie.lines!(mid, range_sweeps, E_pot_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls) *(Ls^2) ; 
                linewidth=1.5, color = :black) 
            CairoMakie.lines!(bottom, range_sweeps, E_pot_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls) *(Ls^2) ; 
                linewidth=1.5, color = :black) 
        end
        count+=2;
    end

    axislegend(mid, position=(1, 0))
    axislegend(bottom, position=(1, 0))




    display(fig)

    if save_bool_fig
        CairoMakie.save(joinpath(p, "Therm_E_pot_t0_L_$(gdf[1][1, :L])_b_" * to_string(1/gdf[1][1, :T]) *
        "_U_" * to_string(gdf[1][1, :U]) * ".png"), fig)
    end
end


################
## spin structure factor S_spin
################
fig = Figure(resolution = (800, 800))

Ts=1.0;Ls=4;Us=0.8;
ana=Sspin_A1_XX_analytic(Us, 1/Ts, 0.5*(μ0-(-μ0)), Ls) ;

mid=Axis(fig[1, 1], title=mytitle,   ylabel=L"$S_{\mathrm{spin}}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],   ylabel=L"$S^{ϕ}_{\mathrm{spin}}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
truebottom=Axis(fig[3, 1],  xlabel=L"$sweep*10$", ylabel=L"$S^{(2),ϕ}_{\mathrm{spin}}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1,]

    count=1;
    for n in 1:length(gdf)
        Ts=gdf[n][1, :T]
        Us=gdf[n][1, :U]
        Ls=gdf[n][1, :L]
        δτ=gdf[n][1, :δτ]
    

        CairoMakie.scatterlines!(mid, gdf[n][!,:_th_meas], gdf[n][!,Symbol("S_spin")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)

        CairoMakie.scatterlines!(bottom, gdf_ϕ[n][!,:_th_meas], gdf_ϕ[n][!,Symbol("S_spin_ϕ")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)

        CairoMakie.scatterlines!(truebottom, gdf_ϕ[n][!,:_th_meas], gdf_ϕ[n][!,Symbol("S2_spin_ϕ")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)
        if count==1
            range_sweeps=1:200;
            range_sweeps=range(1, 200, 30)
            CairoMakie.lines!(mid, range_sweeps, Sspin_A1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ; 
                linewidth=1.5, color = :black) 
            CairoMakie.lines!(bottom, range_sweeps, Sspin_A1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls)  ; 
                linewidth=1.5, color = :black) 
            CairoMakie.lines!(truebottom, range_sweeps, Sspin2_A1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ; 
                linewidth=1.5, color = :black) 
        end
        count+=2;
    end
end

axislegend(truebottom, position=(1, 1))




display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Therm_S_spin_t0_L_$(gdf[1][1, :L])_b_" * to_string(1/gdf[1][1, :T]) *
    "_U_" * to_string(gdf[1][1, :U]) * ".png"), fig)
end

################
## spin susceptibility χ_spin
################
fig = Figure(resolution = (800, 800))



mid=Axis(fig[1, 1], title=mytitle,   ylabel=L"$χ_{\mathrm{spin}}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$sweep*10$", ylabel=L"$χ^{ϕ}_{\mathrm{spin}}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1,]

    count=1;
    for n in 1:length(gdf)
        Ts=gdf[n][1, :T]
        Us=gdf[n][1, :U]
        Ls=gdf[n][1, :L]
        δτ=gdf[n][1, :δτ]
    

        CairoMakie.scatterlines!(mid, gdf[n][!,:_th_meas], gdf[n][!,Symbol("χ_spin")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)

        CairoMakie.scatterlines!(bottom, gdf_ϕ[n][!,:_th_meas], gdf_ϕ[n][!,Symbol("χ_spin_ϕ")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)
        if count==1
            range_sweeps=1:200;
            range_sweeps=range(1, 200, 30)
            CairoMakie.lines!(mid, range_sweeps, χspin_A1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls)  ; 
                linewidth=1.5, color = :black) 
            CairoMakie.lines!(bottom, range_sweeps, χspin_A1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ; 
                linewidth=1.5, color = :black) 
        end
        count+=2;
    end
end
axislegend(mid, position=(1, 1))
axislegend(bottom, position=(1, 1))


display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Therm_chi_spin_t0_L_$(gdf[1][1, :L])_b_" * to_string(1/gdf[1][1, :T]) *
    "_U_" * to_string(gdf[1][1, :U]) * ".png"), fig)
end



################
## B1 bilinear  structure factor S_bil_B1
################
fig = Figure(resolution = (800, 800))

mid=Axis(fig[1, 1], title=mytitle,   ylabel=L"$S_{\mathrm{bilinear}}^{B_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],   ylabel=L"$S_{\mathrm{bilinear}}^{ϕ, B_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
truebottom=Axis(fig[3, 1],  xlabel=L"$sweep*10$", ylabel=L"$S_{\mathrm{bilinear}}^{(2),ϕ, B_1}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1,]

    count=1;
    for n in 1:length(gdf)
        Ts=gdf[n][1, :T]
        Us=gdf[n][1, :U]
        Ls=gdf[n][1, :L]
        δτ=gdf[n][1, :δτ]

        CairoMakie.scatterlines!(mid, gdf[n][!,:_th_meas], gdf[n][!,Symbol("S_bil_B1")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)

        CairoMakie.scatterlines!(bottom,  gdf_ϕ[n][!,:_th_meas],  gdf_ϕ[n][!,Symbol("S_bil_B1_ϕ")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)

        CairoMakie.scatterlines!(truebottom,  gdf_ϕ[n][!,:_th_meas],  gdf_ϕ[n][!,Symbol("S2_bil_B1_ϕ")] ;  
            marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
            label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)
        if count==1
            range_sweeps=1:200;
            range_sweeps=range(1, 200, 30)
            CairoMakie.lines!(mid, range_sweeps, Sbil_B1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ; 
                linewidth=1.5, color = :black) 
            CairoMakie.lines!(bottom, range_sweeps, Sbil_B1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls)  ; 
                linewidth=1.5, color = :black) 
            CairoMakie.lines!(truebottom, range_sweeps, Sbil2_B1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls) ; 
                linewidth=1.5, color = :black) 
        end
        count+=2;
    end
end
axislegend(truebottom, position=(1, 0))


display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Therm_S_bil_B1_t0_L_$(gdf[1][1, :L])_b_" * to_string(1/gdf[1][1, :T]) *
    "_U_" * to_string(gdf[1][1, :U]) * ".png"), fig)
end


################
## B1 bilinear  susceptibility χ_bil_B1
################
fig = Figure(resolution = (800, 800))


mid=Axis(fig[1, 1], title=mytitle,   ylabel=L"$χ_{\mathrm{bilinear}}^{B_1}$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
bottom=Axis(fig[2, 1],  xlabel=L"$sweep*10$", ylabel=L"$χ_{\mathrm{bilinear}}^{ϕ, B_1}$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

for dummy in [1,]

    count=1;
    for n in 1:length(gdf)
        Ts=gdf[n][1, :T]
        Us=gdf[n][1, :U]
        Ls=gdf[n][1, :L]
        δτ=gdf[n][1, :δτ]
    

    CairoMakie.scatterlines!(mid, gdf[n][!,:_th_meas], gdf[n][!,Symbol("χ_bil_B1")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
        label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)

    CairoMakie.scatterlines!(bottom, gdf_ϕ[n][!,:_th_meas], gdf_ϕ[n][!,Symbol("χ_bil_B1_ϕ")] ;  
        marker = :xcross, color = colorschemes[:tab20][count], linewidth=lw,
        label = L"$δτ=%$(δτ), β=%$(1/Ts), L=%$(Ls) $", markersize=ms)
        if count==1
            range_sweeps=1:200;
            range_sweeps=range(1, 200, 30)
            CairoMakie.lines!(mid, range_sweeps, χbil_B1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls)  ; 
                linewidth=1.5, color = :black) 
            CairoMakie.lines!(bottom, range_sweeps, χbil_B1_XX_analytic.([Us for _ in range_sweeps], 1/Ts, 0.5*(μ0-(-μ0)), Ls)  ; 
                linewidth=1.5, color = :black) 
        end
        count+=2;
    end
end
axislegend(mid, position=(1, 1))
axislegend(bottom, position=(1, 1))

display(fig)
if save_bool_fig
    CairoMakie.save(joinpath(p, "Therm_chi_bil_B1_t0_L_$(gdf[1][1, :L])_b_" * to_string(1/gdf[1][1, :T]) *
    "_U_" * to_string(gdf[1][1, :U]) * ".png"), fig)
end
