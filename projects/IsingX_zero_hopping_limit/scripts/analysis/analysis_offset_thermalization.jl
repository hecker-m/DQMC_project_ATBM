# include("../../../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, Plots, LinearAlgebra, CSV
include("../../../../src/Analysis_Fcns/load_fcns.jl")
include("../../../../src/Analysis_Fcns/statistics.jl")
include("/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/Analysis_Fcns/analysis_fcns_thermalization.jl")


#######
## initialize DataFrame for measurements of
## correlation/susceptibility/energy/occupation number etc.
####### 
my_keys=[(key=:E_pot, name="pot_energy"), (key=:SDC_A1_Mx_x, name="S_spin"),
(key=:SDS_A1_Mx_x, name="χ_spin"), (key=:NemC_X, name="S_bil_B1"), 
(key=:NemS_X, name="χ_bil_B1"), (key=:A1p_dQ_C_X, name="S_bil_A1p"), 
(key=:A1p_dQ_S_X, name="χ_bil_A1p")
]
my_keys_ϕ=[(key=:egal, name="pot_energy_ϕ"),
(key=:egal, name="h4_ϕ"), (key=:egal, name="h4_OnSite_ϕ"),
(key=:egal, name="S_spin_ϕ"), (key=:egal, name="S2_spin_ϕ"),
(key=:egal, name="χ_spin_ϕ"), (key=:egal, name="S_bil_B1_ϕ"),
(key=:egal, name="S2_bil_B1_ϕ"), 
(key=:egal, name="χ_bil_B1_ϕ"), (key=:egal, name="S_bil_A1p_ϕ"),
(key=:egal, name="S2_bil_A1p_ϕ"), 
(key=:egal, name="χ_bil_A1p_ϕ")
]

df_cols=(_th_meas=Int[], L=Int[], T=Float64[], U=Float64[], B=Int[], δτ=Float64[]);
for n in eachindex(my_keys)
    if haskey(my_keys[n], :name)
        name=my_keys[n].name;
    else
        name=string(my_keys[n].key);
        if haskey(my_keys[n], :q)
            q1=Int(round(my_keys[n].q[1]/π)) ==0 ? 0 : π;
            q2=Int(round(my_keys[n].q[2]/π)) ==0 ? 0 : π;
            name= name * "_$(q1)$(q2)"
        end
    end
    global df_cols=Base.setindex(df_cols, Float64[], Symbol(name))
end
df= DataFrame(df_cols)

dfϕ_cols=(_th_meas=Int[], L=Int[], T=Float64[], U=Float64[], B=Int[], δτ=Float64[]);
for n in eachindex(my_keys_ϕ)
    if haskey(my_keys_ϕ[n], :name)
        name=my_keys_ϕ[n].name;
    else
        name=string(my_keys_ϕ[n].key);
        if haskey(my_keys_ϕ[n], :q)
            q1=Int(round(my_keys_ϕ[n].q[1]/π)) ==0 ? 0 : π;
            q2=Int(round(my_keys_ϕ[n].q[2]/π)) ==0 ? 0 : π;
            name= name * "_$(q1)$(q2)"
        end
    end
    global dfϕ_cols=Base.setindex(dfϕ_cols, Float64[], Symbol(name))
end
df_ϕ= DataFrame(dfϕ_cols)


dqmcs = [mc_0p05, mc_0p025]

L0=4;
# Us=[0.8, ]
# βs=[1, ]
# paras=[(L=L0, β=β0, U=U0, Pe=true) for U0 in Us, β0 in βs][:]

#######
# to save the computed dataframe
#######
save_bool=false;
save_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_zero_hopping_limit/figures/dataframes/"
my_ϵ="_eps_0pc"
save_name="Therm_Disc8_L$(L0)" * my_ϵ 

prefix_folder="Disc8_";
prefix_file="Disc8_t0_IsX_d_";


for (_para, mc) in enumerate(dqmcs)
    L=mc.model.l.Ls[1]
    N=L^2
    U=mc.model.U
    peierls=mc.model.peierls
    β=mc.parameters.beta
    T=1.0/β
    δτ=mc.parameters.delta_tau

    therm=mc.parameters.thermalization
    meas_rate=mc.parameters.measure_rate

    th_measures=div(therm, meas_rate);

    for _th_measures in 1:th_measures

        μs=zeros(Float64, length(my_keys));
        for (n, tuple) in enumerate(my_keys)
            μs[n] =convert_observable(mc, tuple, _th_measures)
        end
        vec=[_th_measures, L, T, U, Int(peierls), δτ]
        for i=1:length(my_keys)
            append!(vec, μs[i])
        end
        push!(df, vec)

        Epot, h4, h4_OS, S_spin, S2_spin, χ_spin, S_bil_B1, S2_bil_B1, 
        χ_bil_B1, S_bil_A1p, S2_bil_A1p, χ_bil_A1p=convert_observable_ϕ(mc, _th_measures, Val(:Q0πQπ0_offset))
        vec=[_th_measures, L, T, U, Int(peierls), δτ]
        append!(vec, [Epot, h4, h4_OS, S_spin, S2_spin, χ_spin, S_bil_B1, S2_bil_B1,
            χ_bil_B1, S_bil_A1p, S2_bil_A1p, χ_bil_A1p])
        push!(df_ϕ, vec)

    end

end

#########
## Saving the DataFrames if necessary
#########

# if save_bool
#     CSV.write(joinpath(save_path, "df_" * save_name * ".csv"), df)
#     CSV.write(joinpath(save_path, "df_OP_" * save_name * ".csv"), df_OP)
# end



