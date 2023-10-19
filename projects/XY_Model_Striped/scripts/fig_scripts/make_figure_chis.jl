using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings, CSV
include("../../../../src/Analysis_Fcns/large_U_fcns.jl")

p_input="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/XY_Model_Striped/scripts/fig_scripts/"
# p=p_output
p="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/XY_Model_Striped/figures/"

L_plot=8;
peierls=true;
er_lw=1.2;
βmax=20;
Umin=1.9; Umax=2.4;nU=20;
t=1.0; μ0=1.0*t; μs = [μ0; -μ0]
muM=1/4*(μs[1]-μs[2]);
Nϕ=2;

load_bool=true;
if load_bool
    load_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/XY_Model_Striped/figures/dataframes/"
    load_names=["L8_eps_0pc",]
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