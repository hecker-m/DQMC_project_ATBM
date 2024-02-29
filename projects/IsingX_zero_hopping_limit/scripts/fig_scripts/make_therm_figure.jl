using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings, CSV

p_input="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_zero_hopping_limit/scripts/fig_scripts/"
# p=p_output

# fig_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/" * 
#         "IsingX_zero_hopping_limit/figures/Disc8_b_10p0_dtau_0p05/"
fig_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/" * 
        "IsingX_zero_hopping_limit/figures/thermalization/"
if !isdir(fig_path)
    mkdir(fig_path)
end
        
L_plot=4;
ϵ_plot=0.0;
μ0=2.0;
T0=1.0;
U0=0.8;
p=fig_path
save_bool_fig=true;
mytitle=L"$\mathrm{discrete}\!\text{-}\!8 \;\;\mathrm{field}, \quad t=0, \quad μ=%$(μ0), \quad T=%$(T0),
\quad U=%$(U0), \quad L=%$(L_plot)$"





peierls=true;
er_lw=1.2;
βmax=20;


 include(p_input * "figs_offset_therm.jl")
 
 
 
