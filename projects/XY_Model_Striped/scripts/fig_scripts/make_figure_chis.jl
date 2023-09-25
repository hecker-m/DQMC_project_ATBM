using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings, CSV

p_input="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped/scripts/fig_scripts/"
# p=p_output
p="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped/figures/"

load_bool=true;
if load_bool
    load_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped/figures/dataframes/"
    load_names=["L8_Eps_0",]
    for (ii, load_name) in enumerate(load_names)
        if ii ==1
        global df= CSV.read(joinpath(load_path, "df_" * load_name * ".csv"), DataFrame)
        global df_OP= CSV.read(joinpath(load_path, "df_OP_" * load_name * ".csv"), DataFrame)
        global df_ϕ= CSV.read(joinpath(load_path, "df_Phi_" * load_name * ".csv"), DataFrame)
        global df_ϕ_OP= CSV.read(joinpath(load_path, "df_Phi_OP_" * load_name * ".csv"), DataFrame)
        global df_Binder= CSV.read(joinpath(load_path, "df_Binder_" * load_name * ".csv"), DataFrame)
        else 
        global df= vcat(df, CSV.read(joinpath(load_path, "df_" * load_name * ".csv"), DataFrame))
        global df_OP= vcat(df_OP, CSV.read(joinpath(load_path, "df_OP_" * load_name * ".csv"), DataFrame))
        global df_ϕ= vcat(df_ϕ, CSV.read(joinpath(load_path, "df_Phi_" * load_name * ".csv"), DataFrame))
        global df_ϕ_OP= vcat(df_ϕ_OP, CSV.read(joinpath(load_path, "df_Phi_OP_" * load_name * ".csv"), DataFrame))
        global df_Binder= vcat(df_Binder, CSV.read(joinpath(load_path, "df_Binder_" * load_name * ".csv"), DataFrame))
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

include(p_input * "fig_chi_charge.jl")

include(p_input * "fig_chi_pair.jl")
#TODO: Sc order parameter

include(p_input * "fig_chi_spin.jl")
#TODO: Magnetic order parameter

 include(p_input * "fig_chi_B1_nematic.jl")

include(p_input * "fig_chi_A1p_bilinear.jl")

include(p_input * "fig_chi_B1p_bilinear.jl")

include(p_input * "fig_magn_Binder.jl")


include(p_input * "fig_tests.jl")