###########################
### SC bilinear ZZyy 'order parameter'
###########################
"""
SC_bil_ZZyy_OP(mc, model; kwargs...)
Measure the bilinear of the superconducting OP Δ^{Z,y} in a FullBinner.
"""
function Δ_Zy_bil_OP(
        mc::DQMC, model::Model; 
        greens_iterator = Greens(),
        lattice_iterator = nothing, 
        flavor_iterator = nothing,
        kernel = pc_spm_wave_kernel,
        kwargs...
    )
    eltype = geltype(mc)
    obs = FullBinner(Vector{Float64})
    temp = Vector{eltype}(undef, 4)
    return DQMCMeasurement(mc, model, greens_iterator, lattice_iterator, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end
function Δ_0y_bil_OP(
    mc::DQMC, model::Model; 
    greens_iterator = Greens(),
    lattice_iterator = nothing, 
    flavor_iterator = nothing,
    kernel = pc_swave_kernel,
    kwargs...
)
eltype = geltype(mc)
obs = FullBinner(Vector{Float64})
temp = Vector{eltype}(undef, 4)
return DQMCMeasurement(mc, model, greens_iterator, lattice_iterator, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end

function Δ_Xy_bil_OP(
    mc::DQMC, model::Model; 
    greens_iterator = Greens(),
    lattice_iterator = nothing, 
    flavor_iterator = nothing,
    kernel = pc_XX_wave_kernel,
    kwargs...
)
eltype = geltype(mc)
obs = FullBinner(Vector{Float64})
temp = Vector{eltype}(undef, 4)
return DQMCMeasurement(mc, model, greens_iterator, lattice_iterator, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end

function Δ_Ysum_bil_OP(
    mc::DQMC, model::Model; 
    greens_iterator = Greens(),
    lattice_iterator = nothing, 
    flavor_iterator = nothing,
    kernel = pc_YYsum_wave_kernel,
    kwargs...
)
eltype = geltype(mc)
obs = FullBinner(Vector{Float64})
temp = Vector{eltype}(undef, 4)
return DQMCMeasurement(mc, model, greens_iterator, lattice_iterator, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end

@bm function measure!(::Nothing, m::DQMCMeasurement{T}, mc::DQMC, 
    packed_greens) where T<:Union{typeof(pc_spm_wave_kernel), 
    typeof(pc_swave_kernel), typeof(pc_YYsum_wave_kernel), typeof(pc_XX_wave_kernel)}
    m.temp .= zero(eltype(m.temp))
    lat=lattice(mc)
    N = length(lat)
    L= lat.Ls[1]
    Bsrctrg2dir = lat[:Bravais_srctrg2dir]::Matrix{Int}
        @inbounds @fastmath for i in eachindex(lat)
            @simd for j in eachindex(lat)
                δ = Bsrctrg2dir[j, i]   #target i=j+δ
                δy, δx=fldmod1(δ, L) .-(1, 1)
                m.temp[1] += m.kernel(mc, mc.model, (i, j), packed_greens, 0)
                m.temp[2] += (-1)^(δx+δy) * m.kernel(mc, mc.model, (i, j), packed_greens, 0)
                m.temp[3] += (-1)^(δx) * m.kernel(mc, mc.model, (i, j), packed_greens, 0)
                m.temp[4] += (-1)^(δy) * m.kernel(mc, mc.model, (i, j), packed_greens, 0)
            end                             
        end
    push!(m.observable, real.(m.temp)/N^2)
    nothing
end




###########################
### Mx_X order parameter
###########################
"""
Mx_X_OP(mc, model; kwargs...)

Measure the magnetic Mx_x order parameter in a FullBinner.
"""
function Mx_X_OP(
        mc::DQMC, model::Model, dir::Symbol; 
        greens_iterator = Greens(),
        lattice_iterator = nothing, wrapper = nothing, 
        flavor_iterator = nothing,
        kernel = if dir == :x; Mx_X_OP_kernel
        elseif dir == :y; Mx_Y_OP_kernel
        elseif dir == :z; Mx_Z_OP_kernel
        else throw(ArgumentError("`dir` must be :x, :y or :z, but is $dir"))
        end,
        kwargs...
    )
    eltype = geltype(mc)
    obs = FullBinner(Vector{Float64})
    temp = Vector{eltype}(undef, 4)
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(mc, model, greens_iterator, li, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end


@inline Base.@propagate_inbounds function Mx_X_OP_kernel(mc, model::TwoBandModel, i, G::_GM{T}, flv) where T
    return Mx_x_kernel(mc, model, i, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function Mx_Y_OP_kernel(mc, model::TwoBandModel, i, G::_GM{T}, flv) where T
    return Mx_y_kernel(mc, model, i, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function Mx_Z_OP_kernel(mc, model::TwoBandModel, i, G::_GM{T}, flv) where T
    return Mx_z_kernel(mc, model, i, G, flv, field(mc))
end

@bm function measure!(::Nothing, m::DQMCMeasurement{T}, mc::DQMC, packed_greens) where T<:Union{typeof(Mx_X_OP_kernel), 
    typeof(Mx_Y_OP_kernel), typeof(Mx_Z_OP_kernel)}
    m.temp .= zero(eltype(m.temp))
    lat=lattice(mc)
    N = length(lat)
    L= lat.Ls[1]
        @inbounds @fastmath for i in eachindex(lat)
            iy, ix=fldmod1(i, L)                                    
            m.temp[1] += m.kernel(mc, mc.model, i, packed_greens, 0)                  # Mx_x_{0}
            m.temp[2] += (-1)^(ix+iy) *m.kernel(mc, mc.model, i, packed_greens, 0)    # Mx_x_{Q1 +Q2}
            m.temp[3] += (-1)^(ix) *m.kernel(mc, mc.model, i, packed_greens, 0)       # Mx_x_Q1
            m.temp[4] += (-1)^(iy) *m.kernel(mc, mc.model, i, packed_greens, 0)       # Mx_x_Q2
        end
    push!(m.observable, real.(m.temp)/N)
    nothing
end


###########################
### charge-X order parameter, ρ^{X0}_Qi
###########################
"""
proxy_A1p_OP(mc, model; kwargs...)

Proxy A1p order parameter in the charge sector.
"""
function CD_X_OP(
        mc::DQMC, model::Model; 
        greens_iterator = Greens(),
        lattice_iterator = nothing, wrapper = nothing, 
        flavor_iterator = nothing,
        kernel = charge_X_OP_kernel,
        kwargs...
    )
    eltype = geltype(mc)
    obs = FullBinner(Vector{Float64})
    temp = Vector{eltype}(undef, 4)
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(mc, model, greens_iterator, li, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end

@inline Base.@propagate_inbounds function charge_X_OP_kernel(i, N, G::_GM{T}, mc) where T
    return charge_X_OP_kernel(i, N, G, field(mc))
end
@inline Base.@propagate_inbounds function charge_X_OP_kernel(i, N, G::_GM{<: Matrix}, ::AbstractField)
    return -G.val[i, i + 2N] - G.val[i + N, i + 3N] - G.val[i + 2N, i] - G.val[i + 3N, i + N]
end
@inline Base.@propagate_inbounds function charge_X_OP_kernel(i, N, G::_GM{<: Matrix}, 
    ::Union{Discrete_MBF1_X, Discrete_MBF1_X_symm})
    return 0
end
@bm function measure!(::Nothing, m::DQMCMeasurement{typeof(charge_X_OP_kernel)}, mc::DQMC, packed_greens)
    m.temp .= zero(eltype(m.temp))
    lat=lattice(mc)
    N = length(lat)
    L= lat.Ls[1]
        @inbounds @fastmath for i in eachindex(lat)
            iy, ix=fldmod1(i, L)                                    
            m.temp[1] += m.kernel(i, N, packed_greens, mc)                  # ρ0_{0}
            m.temp[2] += (-1)^(ix+iy) *m.kernel(i, N, packed_greens, mc)    # ρ0_{Q1 +Q2}
            m.temp[3] += (-1)^(ix) *m.kernel(i, N, packed_greens, mc)       # ρ0_Q1
            m.temp[4] += (-1)^(iy) *m.kernel(i, N, packed_greens, mc)       # ρ0_Q2
        end
    push!(m.observable, imag.(m.temp)/N)
    nothing
end


###########################
### charge proxy A1p order parameter
###########################
"""
proxy_A1p_OP(mc, model; kwargs...)

Proxy A1p order parameter in the charge sector.
"""
function proxy_A1p_OP(
        mc::DQMC, model::Model; 
        greens_iterator = Greens(),
        lattice_iterator = nothing, wrapper = nothing, 
        flavor_iterator = nothing,
        kernel = proxy_A1p_OP_kernel,
        kwargs...
    )
    eltype = geltype(mc)
    obs = FullBinner(Vector{Float64})
    temp = Vector{eltype}(undef, 3)
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(mc, model, greens_iterator, li, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end

@inline Base.@propagate_inbounds function proxy_A1p_OP_kernel(i, N, G::_GM{T}, mc) where T
    return proxy_A1p_OP_kernel(i, N, G, field(mc))
end
@inline Base.@propagate_inbounds function proxy_A1p_OP_kernel(i, N, G::_GM{<: Matrix}, ::AbstractField)
    return 4 - G.val[i, i] - G.val[i + N, i + N] - G.val[i + 2N, i + 2N] - G.val[i + 3N, i + 3N]
end
@inline Base.@propagate_inbounds function proxy_A1p_OP_kernel(i, N, G::_GM{<: Matrix}, ::Discrete_MBF1_X_symm)
    return 4 - 2*real(G.val[i, i] +G.val[i + N, i + N])
end
@bm function measure!(::Nothing, m::DQMCMeasurement{typeof(proxy_A1p_OP_kernel)}, mc::DQMC, packed_greens)
    m.temp .= zero(eltype(m.temp))
    lat=lattice(mc)
    N = length(lat)
    L= lat.Ls[1]
        @inbounds @fastmath for i in eachindex(lat)
            iy, ix=fldmod1(i, L)                                    
            m.temp[1] += (-1)^(ix+iy) *m.kernel(i, N, packed_greens, mc)    # ρ0_{Q1 +Q2}
            m.temp[2] += (-1)^(ix) *m.kernel(i, N, packed_greens, mc)       # ρ0_Q1
            m.temp[3] += (-1)^(iy) *m.kernel(i, N, packed_greens, mc)       # ρ0_Q2
        end
    push!(m.observable, real.(m.temp)/N)
    nothing
end

###########################
### charge proxy B1 order parameter
###########################

"""
proxy_B1_OP(mc, model; kwargs...)

Proxy B1 order parameter in the charge sector.
"""
function proxy_B1_OP(
        mc::DQMC, model::Model; 
        greens_iterator = Greens(),
        lattice_iterator = nothing, wrapper = nothing, 
        flavor_iterator = nothing,
        kernel = proxy_B1_OP_kernel,
        kwargs...
    )
    eltype = geltype(mc)
    obs = FullBinner(Float64)
    temp = Vector{eltype}(undef, 1)
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(mc, model, greens_iterator, li, flavor_iterator, kernel; obs = obs, temp = temp, kwargs...)
end


@inline Base.@propagate_inbounds function proxy_B1_OP_kernel(is::NTuple{5}, N, G::_GM{T}, mc) where T
    return proxy_B1_OP_kernel(is, N, G, field(mc))
end
@inline Base.@propagate_inbounds function proxy_B1_OP_kernel(is::NTuple{5}, N, G::_GM{<: Matrix}, ::AbstractField)
    i, iPa1, iMa1, iPa2, iMa2 =is
    return -G.val[iMa1, i] + G.val[iMa2, i] - G.val[iPa1, i] + G.val[iPa2, i] - G.val[iMa1 + N, i + N] + G.val[iMa2 + N, i + N] - G.val[iPa1 + N, i + N] + 
    G.val[iPa2 + N, i + N] - G.val[iMa1 + 2N, i + 2N] + G.val[iMa2 + 2N, i + 2N] - G.val[iPa1 + 2N, i + 2N] + G.val[iPa2 + 2N, i + 2N] - 
    G.val[iMa1 + 3N, i + 3N] + G.val[iMa2 + 3N, i + 3N] - G.val[iPa1 + 3N, i + 3N] + G.val[iPa2 + 3N, i + 3N]
end
@inline Base.@propagate_inbounds function proxy_B1_OP_kernel(is::NTuple{5}, N, G::_GM{<: Matrix}, ::Discrete_MBF1_X_symm)
    i, iPa1, iMa1, iPa2, iMa2 =is
    return 2 *real(-G.val[iMa1, i] + G.val[iMa2, i] - G.val[iPa1, i] + G.val[iPa2, i] - 
    G.val[iMa1 + N, i + N] + G.val[iMa2 + N, i + N] - 
    G.val[iPa1 + N, i + N] + G.val[iPa2 + N, i + N])
end

@bm function measure!(::Nothing, m::DQMCMeasurement{typeof(proxy_B1_OP_kernel)}, mc::DQMC, packed_greens)
    m.temp .= zero(eltype(m.temp))
    lat=lattice(mc)
    N = length(lat)
    L= lat.Ls[1]
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        @inbounds @fastmath for i in eachindex(lat)
            iPa1=Bsrcdir2trg[i, 2]           # i + a1
            iMa1=Bsrcdir2trg[i, 1 + L-1]     # i - a1
            iPa2=Bsrcdir2trg[i, 1 + L]       # i + a2
            iMa2=Bsrcdir2trg[i,1+L*L-L]      # i - a2

            m.temp[1] += m.kernel((i, iPa1, iMa1, iPa2, iMa2), N, packed_greens, mc)
        end
    push!(m.observable, real(m.temp[1])/2/N)
    nothing
end

###########################
### lattice iterator for B1 order parameter
###########################

struct EachSitePair_B1_OP <: DirectLatticeIterator end
EachSitePair_B1_OP(::MonteCarloFlavor) = EachSitePair_B1_OP()
output_size(::EachSitePair_B1_OP, l::Lattice) = (1, )
_length(::EachSitePair_B1_OP, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_B1_OP, measurement, mc::DQMC, 
    packed_greens, weight = 1.0
    )
    @timeit_debug "apply!(::EachSitePair_B1_OP, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L= lat.Ls[1]
                    @inbounds @fastmath for px in 1:L, py in 1:L
                        if isodd(px+py)
                            Δp =(-1)^(px-1) -(-1)^(py-1)    #To convert p to direction we effectively compute px, py =(px, py) -(1,1)
                            p=px+(py-1)*L

                          @simd for k in eachindex(lat)
                                kPp=Bsrcdir2trg[k, p]           # k + p                                 
                                temp[1] += weight *Δp *measurement.kernel(mc, mc.model, (k, kPp), packed_greens, 4)
                          end
                        end
                    end
    end
    return 
end

@inline function finalize_temp!(::EachSitePair_B1_OP, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_B1_OP, m) 
    push!(m.observable, real(m.temp[1]))
end
###########################
### lattice iterator for A1p order parameter
###########################

struct EachSitePair_A1p_OP <: DirectLatticeIterator end
EachSitePair_A1p_OP(::MonteCarloFlavor) = EachSitePair_A1p_OP()
output_size(::EachSitePair_A1p_OP, l::Lattice) = (1, )
_length(::EachSitePair_A1p_OP, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_A1p_OP, measurement, mc::DQMC, 
    packed_greens, weight = 1.0
    )
    @timeit_debug "apply!(::EachSitePair_A1p_OP, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L= lat.Ls[1]
                    @inbounds @fastmath for px in 1:L, py in 1:L
                        if iseven(px+py)
                            Δp =(-1)^(px-1) +(-1)^(py-1)    #To convert p to direction we effectively compute px, py =(px, py) -(1,1)
                            p=px+(py-1)*L

                          @simd for k in eachindex(lat)
                                kPp=Bsrcdir2trg[k, p]           # k + p        
                                ky, kx=fldmod1(k, L)                          
                                temp[1] += weight *Δp * (-1)^(kx+ky) * measurement.kernel(mc, mc.model, (k, kPp), packed_greens, 4)
                          end
                        end
                    end
    end
    return 
end

@inline function finalize_temp!(::EachSitePair_A1p_OP, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_A1p_OP, m) 
    push!(m.observable, real(m.temp[1]))
end

###########################
### lattice iterator for B1p order parameter
###########################

struct EachSitePair_B1p_OP <: DirectLatticeIterator end
EachSitePair_B1p_OP(::MonteCarloFlavor) = EachSitePair_B1p_OP()
output_size(::EachSitePair_B1p_OP, l::Lattice) = (1, )
_length(::EachSitePair_B1p_OP, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_B1p_OP, measurement, mc::DQMC, 
    packed_greens, weight = 1.0
    )
    @timeit_debug "apply!(::EachSitePair_B1p_OP, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L= lat.Ls[1]
                    @inbounds @fastmath for px in 1:L, py in 1:L
                        if isodd(px+py)
                            Δp =(-1)^(px-1) -(-1)^(py-1)    #To convert p to direction we effectively compute px, py =(px, py) -(1,1)
                            p=px+(py-1)*L

                          @simd for k in eachindex(lat)
                                kPp=Bsrcdir2trg[k, p]           # k + p        
                                ky, kx=fldmod1(k, L)                          
                                temp[1] += weight *Δp * (-1)^(kx+ky) * measurement.kernel(mc, mc.model, (k, kPp), packed_greens, 4)
                          end
                        end
                    end
    end
    return 
end

@inline function finalize_temp!(::EachSitePair_B1p_OP, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_B1p_OP, m) 
    push!(m.observable, real(m.temp[1]))
end




###########################
### nematic_OP
###########################
function nematic_OP_measurement(
    dqmc::DQMC, model::Model,  greens_iterator; 
        lattice_iterator = EachSitePair_B1_OP(), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = full_nematic_OP_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end

"""
nematic_OP(mc, model; kwargs...)

"""
nematic_OP(args...; kwargs...) = nematic_OP_measurement(args..., Greens(); kwargs...)

###########################
### A1` double-Q order parameter
###########################
function A1p_OP_measurement(
    dqmc::DQMC, model::Model,  greens_iterator; 
        lattice_iterator = EachSitePair_A1p_OP(), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = full_nematic_OP_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end

"""
A1p_OP(mc, model; kwargs...)

"""
A1p_OP(args...; kwargs...) = A1p_OP_measurement(args..., Greens(); kwargs...)

###########################
### B1` double-Q order parameter
###########################

function B1p_OP_measurement(
    dqmc::DQMC, model::Model,   dir::Symbol, greens_iterator; 
        lattice_iterator = EachSitePair_B1p_OP(), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = if dir == :x; full_B1p_x_OP_kernel
        elseif dir == :y; full_B1p_y_OP_kernel
        elseif dir == :z; full_B1p_z_OP_kernel
        else throw(ArgumentError("`dir` must be :x, :y or :z, but is $dir"))
        end,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end

"""
B1p_OP(mc, model; kwargs...)

"""
B1p_OP(args...; kwargs...) = B1p_OP_measurement(args..., Greens(); kwargs...)


###########################
### nematic order parameter kernel
###########################
"""
Calculates the nematic order parameter kernel 
"""
@inline Base.@propagate_inbounds function full_nematic_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_nematic_OP_kernel(mc, model, kkPp, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_nematic_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_nematic_OP_kernel(mc, model, kkPp, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function full_nematic_OP_kernel(
        mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
    )
    k, kPp = kkPp   # k, k+p
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return G0l.val[k + 3N, kPp + 2N]*Gl0.val[kPp, k + N] - (G0l.val[k, kPp + 2N] + 2*G0l.val[k + N, kPp + 3N])*Gl0.val[kPp, k + 2N] + 
    G0l.val[k + N, kPp + 2N]*Gl0.val[kPp, k + 3N] + G0l.val[k + 2N, kPp + 3N]*Gl0.val[kPp + N, k] - 
    G0l.val[k + 3N, kPp + 3N]*(2*Gl0.val[kPp, k] + Gl0.val[kPp + N, k + N]) - G0l.val[k + 2N, kPp + 2N]*(Gl0.val[kPp, k] + 2*Gl0.val[kPp + N, k + N]) + 
    G0l.val[k, kPp + 3N]*Gl0.val[kPp + N, k + 2N] - (2*G0l.val[k, kPp + 2N] + G0l.val[k + N, kPp + 3N])*Gl0.val[kPp + N, k + 3N] - 
    (G0l.val[k + 2N, kPp] + 2*G0l.val[k + 3N, kPp + N])*Gl0.val[kPp + 2N, k] + G0l.val[k + 3N, kPp]*Gl0.val[kPp + 2N, k + N] - 
    (G0l.val[k, kPp] + 2*G0l.val[k + N, kPp + N])*Gl0.val[kPp + 2N, k + 2N] + G0l.val[k + N, kPp]*Gl0.val[kPp + 2N, k + 3N] + 
    G0l.val[k + 2N, kPp + N]*Gl0.val[kPp + 3N, k] - (2*G0l.val[k + 2N, kPp] + G0l.val[k + 3N, kPp + N])*Gl0.val[kPp + 3N, k + N] + 
    G0l.val[k, kPp + N]*Gl0.val[kPp + 3N, k + 2N] - (2*G0l.val[k, kPp] + G0l.val[k + N, kPp + N])*Gl0.val[kPp + 3N, k + 3N] + 
    2*(G00.val[k + N, k + 2N] + G00.val[k + 3N, k])*(Gll.val[kPp, kPp + 3N] + Gll.val[kPp + 2N, kPp + N]) + 
    2*(G00.val[k, k + 3N] + G00.val[k + 2N, k + N])*(Gll.val[kPp + N, kPp + 2N] + Gll.val[kPp + 3N, kPp]) + 
    (G00.val[k, k + 2N] - G00.val[k + N, k + 3N] + G00.val[k + 2N, k] - G00.val[k + 3N, k + N])*(Gll.val[kPp, kPp + 2N] - Gll.val[kPp + N, kPp + 3N] + 
      Gll.val[kPp + 2N, kPp] - Gll.val[kPp + 3N, kPp + N]) + 3*id*(Gl0.val[kPp, k] + Gl0.val[kPp + N, k + N] + Gl0.val[kPp + 2N, k + 2N] + 
      Gl0.val[kPp + 3N, k + 3N])*I[k, kPp]
end
@inline Base.@propagate_inbounds function full_nematic_OP_kernel(
    mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_X_symm
)
    k, kPp = kkPp   # k, k+p
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return 4*real((conj(G00.val[k + N, k]) + G00.val[k, k + N])*(conj(Gll.val[kPp, kPp + N]) + Gll.val[kPp + N, kPp]))+
        6*id*I[k, kPp]*real(Gl0.val[kPp, k] + Gl0.val[kPp + N, k + N])+
        real(2*(conj(G0l.val[k, kPp + N])*Gl0.val[kPp, k + N] + conj(G0l.val[k + N, kPp])*Gl0.val[kPp + N, k] - Gl0.val[kPp + N, k + N]*(G0l.val[k, kPp] + 
        2*real(G0l.val[k, kPp])) - Gl0.val[kPp, k]*(G0l.val[k + N, kPp + N] + 2*real(G0l.val[k + N, kPp + N]))))
end
###########################
### B1` x-order parameter kernel
###########################
"""
Calculates the B1` order parameter kernel in x direction 
"""
@inline Base.@propagate_inbounds function full_B1p_x_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_B1p_x_OP_kernel(mc, model, kkPp, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_B1p_x_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_B1p_x_OP_kernel(mc, model, kkPp, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function full_B1p_x_OP_kernel(
        mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
    )
    k, kPp = kkPp   # k, k+p
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return (-1im)*(-(G0l.val[k + 2N, kPp + 2N]*(Gl0.val[kPp, k + N] + Gl0.val[kPp + N, k])) - G0l.val[k + 3N, kPp + 3N]*(Gl0.val[kPp, k + N] + Gl0.val[kPp + N, k]) + 
    (G0l.val[k + 2N, kPp + 3N] + G0l.val[k + 3N, kPp + 2N])*(Gl0.val[kPp, k] + Gl0.val[kPp + N, k + N]) - 
    (G0l.val[k, kPp + 2N] + G0l.val[k + N, kPp + 3N])*(Gl0.val[kPp, k + 3N] + Gl0.val[kPp + N, k + 2N]) + 
    (G0l.val[k, kPp + 3N] + G0l.val[k + N, kPp + 2N])*(Gl0.val[kPp, k + 2N] + Gl0.val[kPp + N, k + 3N]) - 
    (G0l.val[k + 2N, kPp] + G0l.val[k + 3N, kPp + N])*(Gl0.val[kPp + 2N, k + N] + Gl0.val[kPp + 3N, k]) + 
    (G0l.val[k + 2N, kPp + N] + G0l.val[k + 3N, kPp])*(Gl0.val[kPp + 2N, k] + Gl0.val[kPp + 3N, k + N]) + 
    (G0l.val[k, kPp + N] + G0l.val[k + N, kPp])*(Gl0.val[kPp + 2N, k + 2N] + Gl0.val[kPp + 3N, k + 3N]) - 
    (G00.val[k, k + 2N] - G00.val[k + N, k + 3N] + G00.val[k + 2N, k] - G00.val[k + 3N, k + N])*(Gll.val[kPp, kPp + 3N] - Gll.val[kPp + N, kPp + 2N] + 
      Gll.val[kPp + 2N, kPp + N] - Gll.val[kPp + 3N, kPp]) + (G00.val[k, k + 3N] - G00.val[k + N, k + 2N] + G00.val[k + 2N, k + N] - G00.val[k + 3N, k])*
     (Gll.val[kPp, kPp + 2N] - Gll.val[kPp + N, kPp + 3N] + Gll.val[kPp + 2N, kPp] - Gll.val[kPp + 3N, kPp + N]) + 
    2*id*(Gl0.val[kPp, k + N] + Gl0.val[kPp + N, k])*I[k, kPp] - (Gl0.val[kPp + 2N, k + 3N] + Gl0.val[kPp + 3N, k + 2N])*
     (G0l.val[k, kPp] + G0l.val[k + N, kPp + N] - 2*id*I[k, kPp]))
end

@inline Base.@propagate_inbounds function full_B1p_x_OP_kernel(
    mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, 
    ::Union{Discrete_MBF1_X_symm, Discrete_MBF1_X}
    )
    k, kPp = kkPp   # k, k+p
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return 0
end

###########################
### B1` y-order parameter kernel
###########################
"""
Calculates the B1` order parameter kernel in y direction 
"""
@inline Base.@propagate_inbounds function full_B1p_y_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_B1p_y_OP_kernel(mc, model, kkPp, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_B1p_y_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_B1p_y_OP_kernel(mc, model, kkPp, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function full_B1p_y_OP_kernel(
        mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
    )
    k, kPp = kkPp   # k, k+p
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return  G0l.val[k + 2N, kPp + 2N]*(-Gl0.val[kPp, k + N] + Gl0.val[kPp + N, k]) + G0l.val[k + 3N, kPp + 3N]*(-Gl0.val[kPp, k + N] + Gl0.val[kPp + N, k]) + 
    (G0l.val[k + 2N, kPp + 3N] - G0l.val[k + 3N, kPp + 2N])*(Gl0.val[kPp, k] + Gl0.val[kPp + N, k + N]) - (G0l.val[k, kPp + 2N] + G0l.val[k + N, kPp + 3N])*
     (Gl0.val[kPp, k + 3N] - Gl0.val[kPp + N, k + 2N]) + (G0l.val[k, kPp + 3N] - G0l.val[k + N, kPp + 2N])*(Gl0.val[kPp, k + 2N] + Gl0.val[kPp + N, k + 3N]) - 
    (G0l.val[k + 2N, kPp] + G0l.val[k + 3N, kPp + N])*(Gl0.val[kPp + 2N, k + N] - Gl0.val[kPp + 3N, k]) + (G0l.val[k + 2N, kPp + N] - G0l.val[k + 3N, kPp])*
     (Gl0.val[kPp + 2N, k] + Gl0.val[kPp + 3N, k + N]) + (G0l.val[k, kPp + N] - G0l.val[k + N, kPp])*(Gl0.val[kPp + 2N, k + 2N] + Gl0.val[kPp + 3N, k + 3N]) - 
    (G00.val[k, k + 2N] - G00.val[k + N, k + 3N] + G00.val[k + 2N, k] - G00.val[k + 3N, k + N])*(Gll.val[kPp, kPp + 3N] + Gll.val[kPp + N, kPp + 2N] + Gll.val[kPp + 2N, kPp + N] + 
      Gll.val[kPp + 3N, kPp]) + (G00.val[k, k + 3N] + G00.val[k + N, k + 2N] + G00.val[k + 2N, k + N] + G00.val[k + 3N, k])*
     (Gll.val[kPp, kPp + 2N] - Gll.val[kPp + N, kPp + 3N] + Gll.val[kPp + 2N, kPp] - Gll.val[kPp + 3N, kPp + N]) + 2*id*(Gl0.val[kPp, k + N] - Gl0.val[kPp + N, k])*I[k, kPp] - 
    (Gl0.val[kPp + 2N, k + 3N] - Gl0.val[kPp + 3N, k + 2N])*(G0l.val[k, kPp] + G0l.val[k + N, kPp + N] - 2*id*I[k, kPp])
end

@inline Base.@propagate_inbounds function full_B1p_y_OP_kernel(
    mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, 
    ::Union{Discrete_MBF1_X_symm, Discrete_MBF1_X}
    )
    k, kPp = kkPp   # k, k+p
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return 0
end

###########################
### B1` z-order parameter kernel
###########################
"""
Calculates the B1` order parameter kernel in z direction 
"""
@inline Base.@propagate_inbounds function full_B1p_z_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_B1p_z_OP_kernel(mc, model, kkPp, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_B1p_z_OP_kernel(mc::DQMC, model::TwoBandModel, kkPp::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_B1p_z_OP_kernel(mc, model, kkPp, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function full_B1p_z_OP_kernel(
        mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
    )
    k, kPp = kkPp   # k, k+p
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return  (2*1im)*(G0l.val[k + 3N, kPp + 3N]*Gl0.val[kPp, k] + G0l.val[k + N, kPp + 3N]*Gl0.val[kPp, k + 2N] - G0l.val[k + 2N, kPp + 2N]*Gl0.val[kPp + N, k + N] - 
    G0l.val[k, kPp + 2N]*Gl0.val[kPp + N, k + 3N] + G0l.val[k + 3N, kPp + N]*Gl0.val[kPp + 2N, k] + G0l.val[k + N, kPp + N]*Gl0.val[kPp + 2N, k + 2N] - 
    G0l.val[k + 2N, kPp]*Gl0.val[kPp + 3N, k + N] - G0l.val[k, kPp]*Gl0.val[kPp + 3N, k + 3N] - (G00.val[k + N, k + 2N] + G00.val[k + 3N, k])*
     (Gll.val[kPp, kPp + 3N] + Gll.val[kPp + 2N, kPp + N]) + (G00.val[k, k + 3N] + G00.val[k + 2N, k + N])*(Gll.val[kPp + N, kPp + 2N] + Gll.val[kPp + 3N, kPp]) + 
    id*(-Gl0.val[kPp, k] + Gl0.val[kPp + N, k + N] - Gl0.val[kPp + 2N, k + 2N] + Gl0.val[kPp + 3N, k + 3N])*I[k, kPp])
end

@inline Base.@propagate_inbounds function full_B1p_z_OP_kernel(
    mc, model::TwoBandModel, kkPp::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv,  ::Discrete_MBF1_X_symm
    )
    k, kPp = kkPp   # k, k+p
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return (-2*1im)*(conj(G0l.val[k + N, kPp + N])*conj(Gl0.val[kPp, k]) - conj(G0l.val[k, kPp])*conj(Gl0.val[kPp + N, k + N]) - G0l.val[k + N, kPp + N]*Gl0.val[kPp, k] + 
    G0l.val[k, kPp]*Gl0.val[kPp + N, k + N] + (conj(G00.val[k, k + N]) + G00.val[k + N, k])*(conj(Gll.val[kPp + N, kPp]) + Gll.val[kPp, kPp + N]) - 
    (conj(G00.val[k + N, k]) + G00.val[k, k + N])*(conj(Gll.val[kPp, kPp + N]) + Gll.val[kPp + N, kPp]) + (2*1im)*id*I[k, kPp]*imag(Gl0.val[kPp, k] - Gl0.val[kPp + N, k + N]))
end

