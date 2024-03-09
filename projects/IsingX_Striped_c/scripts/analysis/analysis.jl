include("../../../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, Plots, LinearAlgebra, CSV
include("../../../../src/Analysis_Fcns/analysis_fcns.jl")
include("../../../../src/Analysis_Fcns/load_fcns.jl")
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
save_name="Disc8_L8" * my_ϵ

prefix_folder="Disc8_";
prefix_file="Disc8_IsX_d_";

#######
## initialize DataFrame for measurements of
## correlation/susceptibility/energy/occupation number etc.
#######
my_keys=[(key=:occ, q=(0,0)), (key=:CDS, q=(0,0)),
    (key=:CDS, q=float.((π,π))), (key=:CDS, q=float.((π,0))),
    (key=:CDS, q=float.((0,π))), 
    (key=:CDSxx, q=(0,0)), (key=:CDSxx, q=float.((π,0))),
    (key=:NemC_X, ), (key=:NemS_X, ), 
    (key=:B1_CDS, q=(0,0)), (key=:B1_CDC, q=(0,0)),
    (key=:A1p_dQ_C_X, ), (key=:A1p_dQ_S_X, ),
    (key=:SDS_A1_Mx_x, ), (key=:SDC_A1_Mx_x, ),
    (key=:B1p_dQ_S, q=(0,0)), (key=:B1p_dQ_C, q=(0,0)),
    (key=:SDS_Mx_x, q=(0,0)), (key=:SDS_Mx_x, q=float.((π,π))),
    (key=:SDS_Mx_x, q=float.((π,0))), (key=:SDS_Mx_x, q=float.((0,π))),
    (key=:SDC_Mx_x, q=float.((π,0))), (key=:SDC_Mx_x, q=float.((0,π))),
    (key=:SDC_Mx_x, q=float.((0,0))),
    (key=:PDS_s, q=float.((0,0))),   (key=:PDS_s, q=float.((π,0))),
    (key=:PDS_spm, q=(0,0)), (key=:PDS_spm, q=float.((π,0))),
    (key=:PDS_XX, q=(0,0)), (key=:PDS_YYzz, q=(0,0)),
    (key=:PDS_XX, q=float.((π,0))), (key=:PDS_YYzz, q=float.((π,0))), 
    (key=:PDS_YYzz, q=float.((0,π))),
    (key=:rho_s, name="ρ_s", idx=4), (key=:rho_s, name="ΛL", idx=1),
    (key=:rho_s, name="ΛT", idx=2), (key=:rho_s, name="Λ0", idx=3),
    (key=:kx, name="k_x"), (key=:occ_summed, name="n"),
    (key=:CDC_summed, name="CDC"),
    (key=(:occ_summed, :CDC_summed), name="compress", fcn=compressibility_fcn, norm=compress_norm),
    (key=:E_pot, name="pot_energy"),
    (key=:E, name="energy"),
    (key=(:E, :heat_cap_h2, :heat_cap_h3, :heat_cap_h4), name="spec_heat", fcn=specific_heat_fcn, norm=spec_heat_norm),
    (key=:heat_cap_h2, name="h2"), (key=:heat_cap_h3, name="h3"), (key=:heat_cap_h4, name="h4"),
    (key=:E_TI, name="energy_TI"),
    (key=(:E_TI, :heat_cap_h2_TI, :heat_cap_h3_TI, :heat_cap_h4_TI), name="spec_heat_TI", fcn=specific_heat_fcn, norm=spec_heat_norm),
    (key=:heat_cap_h2_TI, name="h2_TI"), (key=:heat_cap_h3_TI, name="h3_TI"), 
    (key=:heat_cap_h4_TI, name="h4_TI"),
]

df_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[], δτ=Float64[], μ0=Float64[]);
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
    global df_cols=Base.setindex(df_cols, Float64[], Symbol(name))
    global df_cols=Base.setindex(df_cols, Float64[], Symbol(Δname))
end
df= DataFrame(df_cols)
#######
## initialize DataFrame for
## order parameter measurements
#######
my_OP_keys=[(key=:B1_OP , ), (key=:B1_proxy_OP, ),
    (key=:A1p_OP, ), (key=:A1p_proxy_OP, vector=[1]),
    (key=:A1p_proxy_OP, vector=[2]), (key=:Mx_X_OP, vector=[1]),
    (key=:Mx_X_OP, vector=[2]), (key=:Mx_X_OP, vector=[3, 4]), 
    (key=:Δ_Zy_bil_OP, vector=[1])
]

df_OP_cols=(L=Int[], T=Float64[], U=Float64[], B=Int[], δτ=Float64[], μ0=Float64[]);
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
binner_length=117;

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
        _recorder=false, _th_meas=false, _meas=true);
    end

    n_workers != Nworker ? println("Only $(n_workers) finished for parameters $(_para) !!!!") : nothing ;
    println("loaded set $(_para)")
    #global mc_inst=dqmcs[1]     #For later coding, it helps to have one `mc` instance to work with
    μ0=dqmcs[1].model.μs[1];
    δτ=dqmcs[1].parameters.delta_tau;

    ############
    ## Computing the mean and std_errors of all observables 
    ## specified in my_keys, and pushing them into the DataFrame
    ############
    μs=zeros(Float64, length(my_keys));
    σs=zeros(Float64, length(my_keys));
    for (n, tuple) in enumerate(my_keys)
        μs[n], σs[n] =get_value_observable(dqmcs, tuple, my_lvl, binner_length)
    end
    println("done with set $(_para)a")
    vec=[L, T, U, Int(peierls), δτ, μ0]
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
    vec_OP=[L, T, U, Int(peierls), δτ, μ0]
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
