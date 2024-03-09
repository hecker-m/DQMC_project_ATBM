# include("../../../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
# using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, Plots, LinearAlgebra, CSV
include("../../../../src/Analysis_Fcns/analysis_fcns.jl")
include("../../../../src/Analysis_Fcns/load_fcns.jl")
include("../../../../src/Analysis_Fcns/analysis_fcns_phi.jl")
include("../../../../src/Analysis_Fcns/statistics.jl")



Us=[0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4]
βs=[1, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 10, 20]

paras=[(L=8, β=β0, U=U0, Pe=true) for U0 in Us, β0 in βs][:]

δτ_load=0.05;

#######
# to save the computed dataframe
#######
save_bool=true;
save_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped_c/figures/dataframes/"
my_ϵ="_eps_0pc"
save_name="L8" * my_ϵ

prefix_folder="Disc8_";
prefix_file="Disc8_IsX_d_";
#######
## initialize DataFrame for measurements of
## correlation/susceptibility/energy/occupation number etc.
#######
my_keys=[(key=:SDC_A1, q=(0,0)), (key=:SDS_A1, q=(0,0)), 
    (key=:NemC, ), (key=:NemS, ), 
    (key=:A1p_dQ_C, ), (key=:A1p_dQ_S, ), 
    (key=:h4ϕ, ),
    (key=:E_potϕ,),(key=:Sspin_2_A1_xx, ), 
    (key=:Sbil_2_B1_xx, ),(key=:Sbil_2_A1p_xx, )
]



df_ϕ_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[], δτ=Float64[], μ0=Float64[]);
for n in eachindex(my_keys)
    if haskey(my_keys[n], :name)
        name=my_keys[n].name;
        Δname="Δ" *my_keys[n].name ;
    else
        name=string(my_keys[n].key);
        Δname="Δ" *string(my_keys[n].key) ;
        if haskey(my_keys[n], :q)
            q1=Int(round(my_keys[n].q[1]/π)) ==0 ? 0 : π;
            q2=Int(round(my_keys[n].q[2]/π)) ==0 ? 0 : π;
            name= name * "_$(q1)$(q2)"
            Δname= Δname * "_$(q1)$(q2)"
        end
    end
    global df_ϕ_cols=Base.setindex(df_ϕ_cols, Float64[], Symbol(name))
    global df_ϕ_cols=Base.setindex(df_ϕ_cols, Float64[], Symbol(Δname))
end
df_ϕ= DataFrame(df_ϕ_cols)
#######
## initialize DataFrame for
## order parameter measurements
#######
my_OP_keys=[(key=:Mx_OP , ), (key=:B1_OP , ),
    (key=:A1p_OP, )
]

df_ϕ_OP_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[], δτ=Float64[], μ0=Float64[]);
for n in eachindex(my_OP_keys)
    if haskey(my_OP_keys[n], :vector)
        suff="";
        for i in 1:length(my_OP_keys[n].vector)
            suff=suff * "$(my_OP_keys[n].vector[i])"
        end
    else
        suff="";
    end
    name=string(my_OP_keys[n].key) *suff;
    Δname="Δ" *string(my_OP_keys[n].key) *suff;
    global df_ϕ_OP_cols=Base.setindex(df_ϕ_OP_cols, Float64[], Symbol(name))
    global df_ϕ_OP_cols=Base.setindex(df_ϕ_OP_cols, Float64[], Symbol(Δname))
end
df_ϕ_OP= DataFrame(df_ϕ_OP_cols)

#######
## initialize DataFrame for
## Binder cumulant measurements
#######
Binder_keys=[(key=:Binder_magnetic , numbers=(2, 4) ), (key=:Binder_nematic , numbers=(6, 8) ),
(key=:Binder_A1p , numbers=(10, 12) )]

if length(Binder_keys)>0
    df_Binder_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[], δτ=Float64[], μ0=Float64[]);
    for n in eachindex(Binder_keys)
        name=string(Binder_keys[n].key) ;
        Δname="Δ" *string(Binder_keys[n].key);
        global df_Binder_cols=Base.setindex(df_Binder_cols, Float64[], Symbol(name))
        global df_Binder_cols=Base.setindex(df_Binder_cols, Float64[], Symbol(Δname))
    end
    df_Binder= DataFrame(df_Binder_cols)
end

#######
## load all data and evaluate it
#######
my_lvl=7;
binner_length=117
tot_sweeps=2^13;rec_rate=10;
println("bin length $(binner_length) gives $(div(tot_sweeps, rec_rate)/binner_length) bins per walker.")
therm = 1000
for _para in eachindex(paras)
    L=paras[_para].L;
    β=paras[_para].β;  T=1/β;
    U=paras[_para].U;
    peierls=paras[_para].Pe;
    haskey(paras[_para], :sweeps) ? sweeps=paras[_para].sweeps : sweeps = 2^(13);
    haskey(paras[_para], :id) ? jobid=paras[_para].id : jobid = 0;
    N=L^2;
    Nworker=10;
    path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped_c/run_saves/";
    dqmcs = []

    @time begin n_workers=_load(dqmcs, L, T, β, U, peierls, therm, sweeps, Nworker, 
        haskey(paras[_para], :β), jobid;
        path=path, prefix_folder=prefix_folder, prefix_file=prefix_file, eps= my_ϵ, δτ="_dt_" * to_string(δτ_load),
        _recorder=true, _th_meas=false, _meas=true);
    end

    n_workers != Nworker ? println("Only $(n_workers) finished for parameters $(_para) !!!!") : nothing ;
    println("loaded set $(_para)")

    μ0=dqmcs[1].model.μs[1];
    δτ=dqmcs[1].parameters.delta_tau;
    ############
    ## Computing the mean and std_errors of all observables 
    ## specified in my_keys, and pushing them into the DataFrame
    ############

    global my_values=get_observable_using_ϕ(dqmcs,  Val(:Q0πQπ0))
    println("done with set $(_para)a")
    vec=[L, T, U, Int(peierls), δτ, μ0]
    for ν in [2, 3, 6, 7, 10, 11, 13, 14, 4, 8, 12]
        append!(vec, [mean(my_values[ν]), std_error(my_values[ν], binner_length)])
    end
    push!(df_ϕ, vec)

    vec_OP=[L, T, U, Int(peierls), δτ, μ0]
    for ν in [1, 5, 9]
        append!(vec_OP, [mean(my_values[ν]), std_error(my_values[ν], binner_length)])
    end
    push!(df_ϕ_OP, vec_OP)

    if length(Binder_keys)>0
        vec_Binder=[L, T, U, Int(peierls), δτ, μ0]
        for n in eachindex(Binder_keys)
            _Binder_key=Binder_keys[n].key
            _S=Binder_keys[n].numbers[1]
            _S2=Binder_keys[n].numbers[2]

            μ_Binder, Δ_JackKnife_Binder = calculate_Binder(my_values[_S].x , my_values[_S2].x , binner_length)
            append!(vec_Binder, [μ_Binder, Δ_JackKnife_Binder])
        end
        push!(df_Binder, vec_Binder)
    end
end


#########
## Saving the DataFrames if necessary
#########

if save_bool
    CSV.write(joinpath(save_path, "df_Phi_" * save_name * ".csv"), df_ϕ)
    CSV.write(joinpath(save_path, "df_Phi_OP_" * save_name * ".csv"), df_ϕ_OP)
    CSV.write(joinpath(save_path, "df_Binder_" * save_name * ".csv"), df_Binder)
end