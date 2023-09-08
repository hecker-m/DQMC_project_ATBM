################
## combining results from different workers
################
function fill_array!(dqmcs, key, μ_array, σ_array, my_lvl::Int, 
    ::Union{Nothing, EachSitePair_B1, EachDoubleSitePairByDistance, EachDoubleSitePairByDistance_Q1Q2, EachDoubleSitePairByDistance_B1p_Q1Q2})
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
function value_observable(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector ; kwargs...)
    μ_val, σ_val = _mean_std_error!(mc, key, μ_vec, σ_vec, mc[key].lattice_iterator ; kwargs...)
    if abs(imag(μ_val))>10^(-4) || abs(imag(σ_val))>10^(-4)
        println("There is something imaginary.")
        @show μ_val, σ_val
    end       
    return real(μ_val), real(σ_val)
end
function _mean_std_error!(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector, 
    ::Union{EachSitePair_B1, EachDoubleSitePairByDistance, EachDoubleSitePairByDistance_Q1Q2, EachDoubleSitePairByDistance_B1p_Q1Q2} ;
     q::NTuple=(0.0, 0.0), norm=nothing)
    N=length(mc.model.l)
    if norm === nothing
        _norm=1.0
    else
        _norm=norm
    end
    return μ_vec[1]*_norm, σ_vec[1]*_norm
end
function _mean_std_error!(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector, ::Nothing ; 
    q::NTuple=(0.0, 0.0), norm=nothing)
    L =mc.model.l.Ls[1]
    N=length(mc.model.l)
    Nvec=length(μ_vec)
    if norm === nothing
        _norm=1.0/N
    else
        _norm=norm
    end
    nflav= Nvec>N ? div(Nvec, N) : 1; #most measurements are already summed over flavors, but e.g. occupation is not
    μ_val=zero(ComplexF64);
    σ_val=zero(ComplexF64);
    for ix=1:L, iy=1:L, flv=1:nflav
        μ_val+=cis(ix*q[1]+iy*q[2])*μ_vec[ix+(iy-1)*L+(flv -1)*N] 
        σ_val+= σ_vec[ix+(iy-1)*L+(flv -1)*N] 
    end
    return μ_val*_norm, σ_val*_norm
end
function _mean_std_error!(mc::DQMC, key::Symbol, μ_vec::Vector, σ_vec::Vector, iter::EachSitePairByDistance ; 
    q::NTuple=(0.0, 0.0), norm=nothing)
    L =mc.model.l.Ls[1]
    N=length(mc.model.l)
    T=1/mc.parameters.beta;
    Nτ=mc.parameters.slices;
    δτ=mc.parameters.delta_tau; #δτ=1/(T*Nτ)
    norm === nothing ? _norm=1.0 : _norm=norm

    Nvec=length(μ_vec)
    μ_val=zero(ComplexF64);σ_val=zero(ComplexF64);

    for k=1:N
        ky, kx=fldmod1(k, L) 
        dir = [kx ,ky] - [1, 1]
        μ_val+=cis(sum(dir .*q))*μ_vec[k] 
        σ_val+= σ_vec[k] 
    end

    # for ix=1:L, iy=1:L
    #     μ_val+=cis((ix-1)*q[1]+(iy-1)*q[2])*μ_vec[ix+(iy-1)*L] 
    #     σ_val+= σ_vec[ix+(iy-1)*L] 
    # end
    return μ_val*_norm, σ_val*_norm
end
################
## main function to extract physical information
################
function get_value_observable(dqmcs, key::Symbol, Nb::Int,  my_lvl::Int;  kwargs...)
    μ_array, σ_array=make_arrays(dqmcs, key,   my_lvl)
    μ_vec, σ_vec=mean_std_error_combined(μ_array, σ_array, Nb)
    μ_val, σ_val=value_observable(dqmcs[1], key, μ_vec, σ_vec ; kwargs...)
    return μ_val, σ_val
end

################
## main function to compute order parameters
################

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


