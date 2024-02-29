###############
### specific observables function
###############
compressibility_fcn(n, nSq) =nSq - n^2
compress_norm(mc) =length(lattice(mc))* mc.parameters.beta

specific_heat_fcn(en, h2, h3, h4) =h2 +h3 +h4 - en^2
spec_heat_norm(mc) =mc.parameters.beta^2


################
## combining results from different workers
################
function fill_array!(dqmcs, key, μ_array, σ_array, my_lvl::Int, 
    ::Union{Nothing, EachSitePair_B1, EachDoubleSitePairByDistance, 
        EachDoubleSitePairByDistance_Q1Q2, EachDoubleSitePairByDistance_B1p_Q1Q2,
        PS_EachBondPairByBravaisDistance, EachWeightedBond, EachSitePair_A1})
    N_worker=size(μ_array)[2]
    for worker =1:N_worker
        μ_array[:,worker] .=mean(dqmcs[worker][key].observable,1);
        σ_array[:,worker] .=std_error(dqmcs[worker][key].observable, my_lvl);
    end
    return nothing
end
function fill_array!(dqmcs, key, μ_array, σ_array, my_lvl::Int, iter::EachSitePairByDistance)
    N_worker=size(μ_array)[2]
    for worker =1:N_worker
        μ_array[:, worker]=mean(dqmcs[worker][key].observable,1)[:,1,1];
        σ_array[:, worker]=std_error(dqmcs[worker][key].observable, my_lvl)[:,1,1];
    end
    return nothing
end

function make_arrays(dqmcs, key::Symbol,  my_lvl::Int) 
    if isa(mean(dqmcs[1][key].observable), Number)
        N_pos=1;
    else
        N_pos=size(mean(dqmcs[1][key].observable))[1];
    end
    N_worker=length(dqmcs)
    G=typeof(first(mean(dqmcs[1][key].observable)))
    μ_array=Array{G}(undef, N_pos, N_worker);
    σ_array=similar(μ_array);
    fill_array!(dqmcs, key, μ_array, σ_array, my_lvl, dqmcs[1][key].lattice_iterator)

    return μ_array, σ_array
end

function mean_std_error_combined(μνs::Array{T}, σνs::Array{T}, Nb::Int) where T
    N_pos, N_worker=size(μνs)
    σ_result=zeros(T, N_pos)
    μ_result=zeros(T, N_pos)

    for i=1:N_pos
        μb=mean(μνs[i,:])
        if N_worker==1
            σ_result[i]=σνs[i,1]
        else
            out=zero(T);
            for ν=1:N_worker
                out+= -μb^2 + μνs[i,ν]^2 +(Nb-1)*σνs[i,ν]^2
            end
            σ_result[i]=sqrt(out/(N_worker*(N_worker*Nb-1)))
        end
        μ_result[i]=μb
    end
    return μ_result, σ_result
end

################
## evaluating the combined arrays for specific observables of interest
################
"""
    value_observable(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector ; kwargs...)

Using the input vector (`μ_vec`, `σ_vec`), it computes the observable corresponding to `key`.\\
It returns real values only; printing errors if the imaginary part is non-zero.\\
The actual evaluation happens in `_mean_std_error!(mc, key, μ_vec, σ_vec, mc[key].lattice_iterator ; kwargs...)`.
"""
function value_observable(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector ; kwargs...)
    μ_val, σ_val = _mean_std_error!(mc, key, μ_vec, σ_vec, mc[key].lattice_iterator ; kwargs...)
    if abs(imag(μ_val))>10^(-4) || abs(imag(σ_val))>10^(-4)
        println("There is something imaginary for $(key).")
        @show μ_val, σ_val
    end       
    return real(μ_val), real(σ_val)
end
function _mean_std_error!(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector, 
    ::Union{EachSitePair_B1, EachDoubleSitePairByDistance, EachDoubleSitePairByDistance_Q1Q2, 
    EachDoubleSitePairByDistance_B1p_Q1Q2, EachWeightedBond, PS_EachBondPairByBravaisDistance, EachSitePair_A1} ; 
    q::NTuple=(0.0, 0.0), norm=1.0, idx=1 ::Int)

    return μ_vec[idx]*norm, σ_vec[idx]*norm
end


function _mean_std_error!(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector, ::Nothing ; 
    q::NTuple=(0.0, 0.0), norm=1.0/length(lattice(mc)))
    L =mc.model.l.Ls[1]
    N=length(lattice(mc))
    Nvec=length(μ_vec)

    nflav= Nvec>N ? div(Nvec, N) : 1; #most measurements are already summed over flavors, but e.g. occupation is not
    μ_val=zero(ComplexF64);
    σ_val=zero(ComplexF64);
    for ix=1:L, iy=1:L, flv=1:nflav
        μ_val+=cis(ix*q[1]+iy*q[2])*μ_vec[ix+(iy-1)*L+(flv -1)*N] 
        σ_val+= σ_vec[ix+(iy-1)*L+(flv -1)*N] 
    end
    return μ_val*norm, σ_val*norm
end
function _mean_std_error!(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector, iter::EachSitePairByDistance ; 
    q::NTuple=(0.0, 0.0), norm=1.0)
    L =mc.model.l.Ls[1]
    N=length(mc.model.l)

    μ_val=zero(ComplexF64);σ_val=zero(ComplexF64);

    for k=1:N
        ky, kx=fldmod1(k, L) 
        dir = [kx ,ky] - [1, 1]
        μ_val+=cis(sum(dir .*q))*μ_vec[k] 
        σ_val+= σ_vec[k] 
    end
    return μ_val*norm, σ_val*norm
end

###############
### auxiliary functions for below
##############
function first_key(keys)
    return isa(keys, Symbol) ? keys : keys[1]
end

function get_norm(mc::DQMC, tuple)
    return haskey(tuple, :norm) ? tuple.norm(mc) : 1.0
end

################
## main function to extract physical information
################

"""
    get_value_observable(dqmcs, key::Symbol, Nb::Int,  my_lvl::Int;  kwargs...)

The computation of an observable corresponding to `key` with binning level `my_lvl` happens in three steps.\\
First, the means (μ) and std_errors (σ) from all walkers are written into an array (`μ_array`, `σ_array`).\\
Second, the array w.r.t. all walkers is collapsed into one vector (`μ_vec`, `σ_vec`), using the rules of error propagation.\\
Third, we evaluate the corresponding observable using `value_observable(dqmcs[1], key, μ_vec, σ_vec ; kwargs...)`.
"""
function get_value_observable(dqmcs, tuple::NamedTuple, my_lvl::Int,
        binner_length::Int ;  kwargs...)
    f_key =first_key(tuple.key)
    binner = dqmcs[1].measurements[f_key].observable
    return get_value_observable(binner, dqmcs, tuple, my_lvl, binner_length;  kwargs...)
end

function get_value_observable(::LogBinner, dqmcs, tuple::NamedTuple, my_lvl::Int,
        binner_length::Int ;  kwargs...)
    f_key =first_key(tuple.key) #It is assumed that LogBinner measurements have a single key only.

    n_meas_pw=length(dqmcs[1][f_key].observable)
    Nb_pw=div(n_meas_pw, 2^my_lvl)  #number of bins per walker
    
    μ_array, σ_array=make_arrays(dqmcs, f_key,   my_lvl)
    μ_vec, σ_vec=mean_std_error_combined(μ_array, σ_array, Nb_pw)


    if haskey(tuple, :q)
        if !haskey(tuple, :idx)
            μ_val, σ_val =value_observable(dqmcs[1], f_key, μ_vec, σ_vec; q=tuple.q, kwargs...)
        else
            μ_val, σ_val =value_observable(dqmcs[1], f_key, μ_vec, σ_vec; q=tuple.q, idx=tuple.idx, kwargs...)
        end
    else
        if !haskey(tuple, :idx)
            μ_val, σ_val =value_observable(dqmcs[1], f_key, μ_vec, σ_vec; kwargs...)
        else
            μ_val, σ_val =value_observable(dqmcs[1], f_key, μ_vec, σ_vec; idx=tuple.idx, kwargs...)
        end
    end
    #μ_val, σ_val=value_observable(dqmcs[1], f_key, μ_vec, σ_vec ; kwargs...)
    return μ_val, σ_val
end

function get_value_observable(binner::FullBinner,dqmcs, tuple::NamedTuple, my_lvl::Int ,
        binner_length::Int ;  kwargs...)

    if isa(tuple.key, Symbol)   
    #for multi-variate quantities (e.g. heat capacity), tuple.key is a tuple of Symbols
        get_value_observable_Scalar(binner ,dqmcs, tuple, my_lvl, binner_length;  kwargs...)
    else
        get_value_observable_JackKnife(binner ,dqmcs, tuple, my_lvl, binner_length;  kwargs...)
    end
end
function get_value_observable_Scalar(binner::FullBinner,dqmcs, tuple::NamedTuple, my_lvl::Int ,
    binner_length::Int ;  norm=get_norm(dqmcs[1], tuple))

    type=eltype(dqmcs[1].measurements[tuple.key].observable.x[1])
    all_values=FullBinner(type)
    for walker in eachindex(dqmcs)
        append!(all_values, dqmcs[walker].measurements[tuple.key].observable.x)
    end

    return mean(all_values) *norm, std_error(all_values, binner_length) *norm
end

function get_value_observable_JackKnife(binner::FullBinner,dqmcs, tuple::NamedTuple, my_lvl::Int ,
        binner_length::Int ;  norm=get_norm(dqmcs[1], tuple))

    keys=tuple.key
    all_vectors = map(keys) do key
        type=typeof(dqmcs[1].measurements[key].observable.x[1])
        full_vector=Vector{type}(undef, 0)
        for walker in eachindex(dqmcs)
            full_vector=vcat(full_vector, dqmcs[walker].measurements[key].observable.x)
        end
        return full_vector
    end
    binned_vectors = prebinning(binner_length, all_vectors...)


    μ_g, σJackKnife_g = jackknife(tuple.fcn, binned_vectors...)


    return μ_g *norm, σJackKnife_g *norm
end
################
## main function to compute order parameters
################

"""
    get_OP_value(dqmcs, key::Symbol, Nb::Int,  my_lvl::Int;  vec=[])

Computes the mean of the absolute value of the measured order parameter (OP).
Note that by default, OPs are saved in FullBinners, such that one can look at histograms.
"""
function get_OP_value(dqmcs, key::Symbol, Nb::Int,  my_lvl::Int;  vec=[])

    N_worker=length(dqmcs)
    G=typeof(first(mean(dqmcs[1][key].observable)))
    obs=FullBinner(G)
    for worker =1:N_worker
        N_meas=length(dqmcs[worker][key].observable)
        for n in 1:N_meas
            OP_push!(obs, dqmcs[worker][key].observable.x[n], vec, Val(length(vec)))
        end
    end
    return mean(obs), std_error(obs, Nb)
end
function OP_push!(obs::FullBinner, obs_val, vec, ::Val{0})
    push!(obs, abs(obs_val))
    return nothing
end

function OP_push!(obs::FullBinner, obs_val, vec, ::Val{1})
    push!(obs, abs(obs_val[vec[1]]))
    return nothing
end
function OP_push!(obs::FullBinner, obs_val, vec, ::Val{2})
    push!(obs, sqrt(abs2(obs_val[vec[1]])+abs2(obs_val[vec[2]])))
    return nothing
end









################
## evaluating observables using HS-fields ϕ
################
#Computes the magnetic susceptibility using the recorded configurations
#of ϕ. It takes quite a while---recommended to not use the entire sample set.
function get_χϕ(dqmcs; N_worker::Int=length(dqmcs), lconfig::Int=length(dqmcs[1].recorder.configs))
    N=length(dqmcs[1].model.l)
    T=1/dqmcs[1].parameters.beta;
    Nτ=dqmcs[1].parameters.slices;
    χϕs=zeros(lconfig*N_worker);
    for worker =1:N_worker
        conf=dqmcs[worker].recorder.configs
        for τ_MC in 1:lconfig
            out=0.0;    
            for ix=1:L, iy=1:L, jx=1:L, jy=1:L, ℓp=1:Nτ, ℓ=1:Nτ
                out+=conf[τ_MC][1,ix+(iy-1)*L,ℓ] * conf[τ_MC][1,jx+(jy-1)*L,ℓp]
            end
            χϕs[τ_MC+(worker-1)*lconfig]=1/(T*N*Nτ^2)*real(out)
        end
    end
    return mean(χϕs), std_error(χϕs)
end

#=
key=:occ
for my_lvl=6:12
    σ=1/N*sum(std_error(dqmcs[1].measurements[key].observable,my_lvl))
    println("bin_size=2^$(my_lvl), σ=", σ)
end

my_lvl=9;
n_meas=length(dqmcs[1][key].observable)
Nb=div(n_meas, 2^my_lvl)



μ_occ, σ_occ=get_value_observable(dqmcs, :occ, Nb,  my_lvl; q=(0, 0))
μ_χSDW, σ_χSDW=get_value_observable(dqmcs, :SDC_Mx_z, Nb,  my_lvl; q=(0, 0))
μ_χSDW_ϕ, σ_χSDW_ϕ=get_χϕ(dqmcs, N_worker=1, lconfig=200)


println("μ_occ = $(μ_occ) ± $(σ_occ)")
println("μ_χSDW = $(μ_χSDW) ± $(σ_χSDW)")
println("Using χSDW(ϕ)=$(μ_χSDW_ϕ), and χSDW=$(μ_χSDW), 
the relation (2U)χSDW+1=χSDW(ϕ) becomes $(μ_χSDW_ϕ)=$((2U)*μ_χSDW+1)")
=#



