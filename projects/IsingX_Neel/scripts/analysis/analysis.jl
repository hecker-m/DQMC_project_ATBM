include("../../../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, Plots, LinearAlgebra, CSV
include("../../../../src/Analysis_Fcns/analysis_fcns.jl")
include("../../../../src/Analysis_Fcns/load_fcns.jl")


# Us=[0.6, 0.7, 0.8, 1.0, 1.2, 1.3]
# βs=[1, 4, 10, 14, 20]
Us=[0.7, ]
βs=[20, ]
paras=[(L=8, β=β0, U=U0, Pe=true) for U0 in Us, β0 in βs][:]

#######
# to save the computed dataframe
#######
save_bool=false;
save_path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Neel/figures/dataframes/"
save_name="L8" 


#######
## initialize DataFrame for measurements of
## correlation/susceptibility/energy/occupation number etc.
#######
my_keys=[(key=:occ, q=(0,0)), (key=:CDS, q=(0,0)),
    (key=:CDS, q=float.((π,π))), 
    (key=:SDS_Mx_x, q=(0,0)), (key=:SDS_Mx_x, q=float.((π,π))),
    (key=:SDS_Mx_y, q=(0,0)), (key=:SDS_Mx_y, q=float.((π,π))),
    (key=:SDS_Mx_z, q=(0,0)), (key=:SDS_Mx_z, q=float.((π,π))),
    (key=:PDS_s, q=float.((0,0))), (key=:PDS_spm, q=(0,0)), 
    (key=:rho_s, q=(0,0), idx=4), (key=:rho_s, q=(0,0), idx=1),
    (key=:rho_s, q=(0,0), idx=2), (key=:rho_s, q=(0,0), idx=3),
    (key=:kx, q=(0,0))]

df_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[]);
for n in eachindex(my_keys)
    name=string(my_keys[n].key);
    Δname="Δ" *string(my_keys[n].key) ;
    if haskey(my_keys[n], :q)
        q1=Int(round(my_keys[n].q[1]/π)) ==0 ? 0 : π;
        q2=Int(round(my_keys[n].q[2]/π)) ==0 ? 0 : π;
        name= name * "_$(q1)$(q2)"
        Δname= Δname * "_$(q1)$(q2)"
    end
    if haskey(my_keys[n], :idx)
        name= name * "_$(my_keys[n].idx)"
        Δname= Δname * "_$(my_keys[n].idx)"
    end
    global df_cols=Base.setindex(df_cols, Float64[], Symbol(name))
    global df_cols=Base.setindex(df_cols, Float64[], Symbol(Δname))
end
df= DataFrame(df_cols)
#######
## initialize DataFrame for
## order parameter measurements
#######
my_OP_keys=[(key=:Mx_X_OP, vector=[1]), (key=:Mx_X_OP, vector=[2]), 
    (key=:Mx_X_OP, vector=[3, 4]), (key=:Mx_Y_OP, vector=[1]), 
    (key=:Mx_Y_OP, vector=[2]),  (key=:Mx_Y_OP, vector=[3, 4]), 
    (key=:Mx_Z_OP, vector=[1]), (key=:Mx_Z_OP, vector=[2]), 
    (key=:Mx_Z_OP, vector=[3, 4])]

df_OP_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[]);
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
    global df_OP_cols=Base.setindex(df_OP_cols, Float64[], Symbol(name))
    global df_OP_cols=Base.setindex(df_OP_cols, Float64[], Symbol(Δname))
end
df_OP= DataFrame(df_OP_cols)


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
    path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Neel/run_saves/";
    global dqmcs = []

    @time begin n_workers=_load(dqmcs, L, T, β, U, peierls, therm, sweeps, Nworker, 
        haskey(paras[_para], :β), jobid;
        path=path, prefix_folder="D2_", prefix_file="DSx_FP_",
        _recorder=false, _th_meas=false, _meas=true);
    end
    n_workers != Nworker ? println("Only $(n_workers) finished for parameters $(_para) !!!!") : nothing ;
    println("loaded set $(_para)")

    ############
    ## Computing the mean and std_errors of all observables 
    ## specified in my_keys, and pushing them into the DataFrame
    ############
    μs=zeros(Float64, length(my_keys));
    σs=zeros(Float64, length(my_keys));
    for (n, tuple) in enumerate(my_keys)
        n_meas=length(dqmcs[1][tuple.key].observable)
        Nb=div(n_meas, 2^my_lvl)
        if !haskey(tuple, :idx)
            μs[n], σs[n] =get_value_observable(dqmcs, tuple.key, Nb,  my_lvl; q=tuple.q)
        else
            μs[n], σs[n] =get_value_observable(dqmcs, tuple.key, Nb,  my_lvl; q=tuple.q, idx=tuple.idx)
        end
        if tuple.key==:CDS && iszero(tuple.q[1]) && iszero(tuple.q[2])
            μs[n]=μs[n]-N/T*μs[n-1]^2 
        end
    end
    println("done with set $(_para)a")
    vec=[L, T, U, Int(peierls)]
    for i=1:length(my_keys)
        append!(vec, [μs[i], σs[i]])
    end
    push!(df, vec)

    ############
    ## Computing the mean and std_errors of all order parameter observables 
    ## specified in my_OP_keys, and pushing them into the DataFrame
    ############
    μs=zeros(Float64, length(my_OP_keys));
    σs=zeros(Float64, length(my_OP_keys));
    for (n, tuple) in enumerate(my_OP_keys)
        n_meas=length(dqmcs[1][tuple.key].observable)
        Nb=div(n_meas, 2^my_lvl)

        haskey(tuple, :vector) ? vector=tuple.vector : vector=[]

        μs[n], σs[n] =get_OP_value(dqmcs, tuple.key, Nb,  my_lvl; vec=vector)
 
    end
    println("done with set $(_para)b")
    vec_OP=[L, T, U, Int(peierls)]
    for i=1:length(my_OP_keys)
        append!(vec_OP, [μs[i], σs[i]])
    end
    push!(df_OP, vec_OP)
end



#########
## Saving the DataFrames if necessary
#########

if save_bool
    CSV.write(joinpath(save_path, "df_" * save_name * ".csv"), df)
    CSV.write(joinpath(save_path, "df_OP_" * save_name * ".csv"), df_OP)
end