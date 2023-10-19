using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings

p_input="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Neel/scripts/fig_scripts/"
# p=p_output
p_screenshot="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Neel/figures/"
p="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Neel/figures/"

df=sort(df, [:U, :T, :B, :L]);
df_LT=sort(unique(df[:,[:T,:L]]),:T, rev=true)
df_LU=sort(unique(df[:,[:U,:L]]),:U)

df_OP=sort(df_OP, [:U, :T, :B, :L]);
df_OP_LT=sort(unique(df_OP[:,[:T,:L]]),:T , rev=true) 
df_OP_LU=sort(unique(df_OP[:,[:U,:L]]),:U) 

# df_ϕ=sort(df_ϕ, [:U, :T, :B, :L]);
# df_ϕ_LT=sort(unique(df_ϕ[:,[:T,:L]]),:T, rev=true)
# df_ϕ_LU=sort(unique(df_ϕ[:,[:U,:L]]),:U)

# df_ϕ_OP=sort(df_ϕ_OP, [:U, :T, :B, :L]);
# df_ϕ_OP_LT=sort(unique(df_ϕ_OP[:,[:T,:L]]),:T , rev=true)
# df_ϕ_OP_LU=sort(unique(df_ϕ_OP[:,[:U,:L]]),:U)


peierls=true;
er_lw=1.2;
βmax=30;

 include(p_input * "fig_chi_charge.jl")

# include(p_input * "fig_chi_pair.jl")

# include(p_input * "fig_chi_spin.jl")

# include(p_input * "fig_magn_Binder.jl")

# include(p_input * "fig_phase_stiffness.jl")

