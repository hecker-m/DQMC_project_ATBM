using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings, CSV

p_input="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_zero_hopping_limit/scripts/fig_scripts/"
# p=p_output

fig_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/" * 
        "IsingX_zero_hopping_limit/figures/Disc8_b_10p0_dtau_0p05/"
if !isdir(fig_path)
    mkdir(fig_path)
end
        
L_plot=8;
ϵ_plot=0.0;
μ0=2.0;
#δτ=0.025;
p=fig_path
save_bool_fig=true;
mytitle=L"$\mathrm{discrete}\!\text{-}\!8 \;\;\mathrm{field}, \quad μ=%$(μ0),\quad t=0,\quad δτ=%$(δτ)$"

load_bool=false;
if load_bool
    load_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/\
            IsingX_Striped_d/figures/dataframes/df_Disc4/df_dt_0p05"
    ϵs=[0.0, 0.0];
    load_names=["L8_eps_0pc", ]

    for (ii, load_name) in enumerate(load_names)
        if ii ==1
            global df= CSV.read(joinpath(load_path, "df_" * load_name * ".csv"), DataFrame)
            df[!,:ϵ0] .=ϵs[ii]
            global df_OP= CSV.read(joinpath(load_path, "df_OP_" * load_name * ".csv"), DataFrame)
            df_OP[!,:ϵ0] .=ϵs[ii]
            global df_ϕ= CSV.read(joinpath(load_path, "df_Phi_" * load_name * ".csv"), DataFrame)
            df_ϕ[!,:ϵ0] .=ϵs[ii]
            global df_ϕ_OP= CSV.read(joinpath(load_path, "df_Phi_OP_" * load_name * ".csv"), DataFrame)
            df_ϕ_OP[!,:ϵ0] .=ϵs[ii]
            global df_Binder= CSV.read(joinpath(load_path, "df_Binder_" * load_name * ".csv"), DataFrame)
            df_Binder[!,:ϵ0] .=ϵs[ii]
        else 
            df_add=CSV.read(joinpath(load_path, "df_" * load_name * ".csv"), DataFrame)
            df_add[!,:ϵ0] .=ϵs[ii]
            global df= vcat(df, df_add)
            df_OP_add=CSV.read(joinpath(load_path, "df_OP_" * load_name * ".csv"), DataFrame)
            df_OP_add[!,:ϵ0] .=ϵs[ii]
            global df_OP= vcat(df_OP, df_OP_add)
            df_ϕ_add=CSV.read(joinpath(load_path, "df_Phi_" * load_name * ".csv"), DataFrame)
            df_ϕ_add[!,:ϵ0] .=ϵs[ii]
            global df_ϕ= vcat(df_ϕ, df_ϕ_add)
            df_ϕ_OP_add=CSV.read(joinpath(load_path, "df_Phi_OP_" * load_name * ".csv"), DataFrame)
            df_ϕ_OP_add[!,:ϵ0] .=ϵs[ii]
            global df_ϕ_OP= vcat(df_ϕ_OP, df_ϕ_OP_add)
            df_Binder_add=CSV.read(joinpath(load_path, "df_Binder_" * load_name * ".csv"), DataFrame)
            df_Binder_add[!,:ϵ0] .=ϵs[ii]
            global df_Binder= vcat(df_Binder, df_Binder_add)
        end
    end
end



#########
## Sorting and Plotting
#########
df=sort(df, [:U, :T, :B, :L]);
df_LT=sort(unique(df[:,[:T,:L]]),:T, rev=true)
df_LU=sort(unique(df[:,[:U,:L]]),:U)

df_OP=sort(df_OP, [:U, :T, :B, :L]);
df_OP_LT=sort(unique(df_OP[:,[:T,:L]]),:T , rev=true) 
df_OP_LU=sort(unique(df_OP[:,[:U,:L]]),:U) 

df_ϕ=sort(df_ϕ, [:U, :T, :B, :L]);
df_ϕ_LT=sort(unique(df_ϕ[:,[:T,:L]]),:T, rev=true)
df_ϕ_LU=sort(unique(df_ϕ[:,[:U,:L]]),:U)

df_ϕ_OP=sort(df_ϕ_OP, [:U, :T, :B, :L]);
df_ϕ_OP_LT=sort(unique(df_ϕ_OP[:,[:T,:L]]),:T , rev=true)
df_ϕ_OP_LU=sort(unique(df_ϕ_OP[:,[:U,:L]]),:U)


peierls=true;
er_lw=1.2;
βmax=20;


 include(p_input * "figs_offset.jl")
 
 
 
