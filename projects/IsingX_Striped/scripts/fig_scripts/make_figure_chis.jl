using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings


p="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/IsingX_Striped/figures/"
df=sort(df, [:U, :T, :B, :L]);
peierls=true;
er_lw=1.2;
βmax=20;

include(p * "scripts/fig_chi_charge.jl")

include(p * "scripts/fig_chi_pair.jl")
#TODO: Sc order parameter

include(p * "scripts/fig_chi_spin.jl")
#TODO: Magnetic order parameter

include(p * "scripts/fig_chi_B1_nematic.jl")

include(p * "scripts/fig_chi_A1p_bilinear.jl")

include(p * "scripts/fig_chi_B1p_bilinear.jl")




