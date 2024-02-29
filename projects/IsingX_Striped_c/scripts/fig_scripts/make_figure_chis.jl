using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings, CSV

path_Analysis_Fcns="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/Analysis_Fcns/"
include(path_Analysis_Fcns * "analytical_fcns.jl")

p_input="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped_c/scripts/fig_scripts/"
# p=p_output

fig_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/" * 
        "IsingX_Striped_c/figures/Disc8/"

L_plot=8;
μ0_plot=1.0;
δτ_plot=0.05;
ϵ_plot=0.0;
if isa(L_plot, Int)
    p=fig_path * "figs_L$(L_plot)/"
    if !iszero(ϵ_plot)
        p=p* "figs_eps_$(Int(100*ϵ_plot))pc/"
    end
else
    p=fig_path
end
if !isdir(p)
    mkdir(p)
end
save_bool_fig=true;

load_bool=false;
if load_bool
    load_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped/figures/dataframes/"
    # load_names=["L8_Eps_0", "L10_Eps_0"]
    ϵs=[0.0,0.05,0.1,0.2];
    load_names=["L8_Eps_0", "L8_eps_5pc", "L8_eps_10pc", "L8_eps_20pc"]

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
tit="";
@isdefined(L_plot) ? tit=tit * ",\\quad L=$(L_plot)" : nothing ;
@isdefined(μ0_plot) ? tit = tit * ",\\quad μ=$(μ0_plot)"  : nothing;
@isdefined(δτ_plot) ? tit = tit * ",\\quad δτ=$(δτ_plot)"  : nothing;

my_title=latexstring("\$\\mathrm{discrete}\\!\\text{-}\\!8 \\;\\;\\mathrm{field}" * tit * "\$")
my_title2=latexstring("\$\\mathrm{discrete}\\!\\text{-}\\!8 \\;\\;\\mathrm{field}:\\quad \\times \\;↔\\;  
\\langle \\hat{c}\\hat{c}.. \\rangle, \\qquad □ \\;↔\\;  \\langle ϕ.. \\rangle" * 
tit * "\$")



minU=findmin(df[!,:U])[1];
maxU=findmax(df[!,:U])[1];
rangeU=range(minU, maxU, 30)

minT=findmin(df[!,:T])[1];
maxT=findmax(df[!,:T])[1];
rangeT=range(minT, maxT, 30)

include(p_input * "fig_chi_charge.jl")

include(p_input * "fig_chi_pair.jl")

include(p_input * "fig_chi_spin.jl")

include(p_input * "fig_chi_B1_nematic.jl")

include(p_input * "fig_chi_A1p_bilinear.jl")

include(p_input * "fig_chi_B1p_bilinear.jl")

# include(p_input * "fig_magn_Binder.jl")
include(p_input * "fig_phase_stiffness.jl")

include(p_input * "fig_spec_heat.jl")

# include(p_input * "fig_tests.jl")
     
# include(p_input * "fig_combined_Binder.jl")

# include(p_input * "fig_strain.jl")