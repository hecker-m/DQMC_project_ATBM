# include("../../../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
# using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, Plots, LinearAlgebra, CSV
include("../../../../src/Analysis_Fcns/analysis_fcns.jl")
include("../../../../src/Analysis_Fcns/load_fcns.jl")
include("../../../../src/Analysis_Fcns/analysis_fcns_phi.jl")


Us=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.4]
βs=[ 1, 2,  6, 10, 20]
paras=[(L=8, β=β0, U=U0, Pe=true) for U0 in Us, β0 in βs][:]


#######
## initialize DataFrame for measurements of
## correlation/susceptibility/energy/occupation number etc.
#######
my_keys=[(key=:SDC_A1, q=(0,0)), (key=:SDS_A1, q=(0,0)), 
    (key=:NemC, q=(0,0)), (key=:NemS, q=(0,0)), 
    (key=:A1p_dQ_C, q=(0,0)), (key=:A1p_dQ_S, q=(0,0)),
    (key=:B1p_dQ_C, q=(0,0)), (key=:B1p_dQ_S, q=(0,0))]

df_ϕ_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[]);
for n in eachindex(my_keys)
    q1=Int(round(my_keys[n].q[1]/π)) ==0 ? 0 : π;
    q2=Int(round(my_keys[n].q[2]/π)) ==0 ? 0 : π;
    name=string(my_keys[n].key) *"_$(q1)$(q2)";
    Δname="Δ" *string(my_keys[n].key) *"_$(q1)$(q2)";
    global df_ϕ_cols=Base.setindex(df_ϕ_cols, Float64[], Symbol(name))
    global df_ϕ_cols=Base.setindex(df_ϕ_cols, Float64[], Symbol(Δname))
end
df_ϕ= DataFrame(df_ϕ_cols)
#######
## initialize DataFrame for
## order parameter measurements
#######
my_OP_keys=[(key=:Mx_OP , ), (key=:B1_OP , ),
    (key=:A1p_OP, ), (key=:B1p_OP, )]

df_ϕ_OP_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[]);
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
## load all data and evaluate it
#######
my_lvl=7;
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
    path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/Heisenberg_Striped/run_saves/";
    dqmcs = []

    @time begin n_workers=_load(dqmcs, L, T, β, U, peierls, therm, sweeps, Nworker, 
        haskey(paras[_para], :β), jobid;
        path=path, prefix_folder="D_", prefix_file="D_Hn_",
        _recorder=true, _th_meas=false, _meas=true);
    
    n_workers != Nworker ? println("Only $(n_workers) finished for parameters $(_para) !!!!") : nothing ;
    println("loaded set $(_para)")

    ############
    ## Computing the mean and std_errors of all observables 
    ## specified in my_keys, and pushing them into the DataFrame
    ############
    my_values=get_observable_using_ϕ(dqmcs,  Val(:Q0πQπ0))
    println("done with set $(_para)a")
    vec=[L, T, U, Int(peierls)]
    for ν in [2, 3, 6, 7, 10, 11, 14, 15]
        append!(vec, [mean(my_values[ν]), std_error(my_values[ν], 2^(my_lvl-1))])
    end
    push!(df_ϕ, vec)

    vec_OP=[L, T, U, Int(peierls)]
    for ν in [1, 5, 9, 13]
        append!(vec_OP, [mean(my_values[ν]), std_error(my_values[ν], 2^(my_lvl-1))])
    end
    push!(df_ϕ_OP, vec_OP)
    end
end

#########
## Saving the DataFrames if necessary
#########
# CSV.write(joinpath(p, "dataframe.csv"), df)
# CSV.write(joinpath(p, "dataframe_OP.csv"), df_OP)