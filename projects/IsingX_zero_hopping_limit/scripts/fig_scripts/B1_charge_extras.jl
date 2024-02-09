

###########################
### Zero-hopping (t=0) B1_charge_density susceptibility kernel
###########################


"""
Calculates the B1_charge_density susceptibility kernel 
"""
@inline Base.@propagate_inbounds function t0_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10}, G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return t0_B1_charge_density_kernel(mc, model, isks, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function t0_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10},  G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return t0_B1_charge_density_kernel(mc, model, isks, (G, G, G, G), flv, field)
end


@inline Base.@propagate_inbounds function t0_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField  )
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 


    return nothing
end

@inline Base.@propagate_inbounds function t0_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::Union{Abstract_DiscreteMBF1_X_symm, Cont_MBF1_X_symm, Discrete_MBF2_symm})
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return (1/N)*real(2*((id - G0l.val[i, i])*Gl0.val[iMa1, iMa1] + (id - G0l.val[i, i])*Gl0.val[iMa2, iMa2] + id*Gl0.val[iPa1, iPa1] - 
      G0l.val[i, i]*Gl0.val[iPa1, iPa1] + id*Gl0.val[iPa2, iPa2] - G0l.val[i, i]*Gl0.val[iPa2, iPa2] - 
      G0l.val[i + N, i]*(Gl0.val[iMa1, iMa1 + N] + Gl0.val[iMa2, iMa2 + N] + Gl0.val[iPa1, iPa1 + N] + Gl0.val[iPa2, iPa2 + N]) + 
      id*Gl0.val[iMa1 + N, iMa1 + N] - G0l.val[i + N, i + N]*Gl0.val[iMa1 + N, iMa1 + N] + id*Gl0.val[iMa2 + N, iMa2 + N] - 
      G0l.val[i + N, i + N]*Gl0.val[iMa2 + N, iMa2 + N] + id*Gl0.val[iPa1 + N, iPa1 + N] - 
      G0l.val[i + N, i + N]*Gl0.val[iPa1 + N, iPa1 + N] - G0l.val[i, i + N]*(Gl0.val[iMa1 + N, iMa1] + Gl0.val[iMa2 + N, iMa2] + 
        Gl0.val[iPa1 + N, iPa1] + Gl0.val[iPa2 + N, iPa2]) + (id - G0l.val[i + N, i + N])*Gl0.val[iPa2 + N, iPa2 + N]))
end
###########################
### Zero-hopping, diagonal-only (t=0) B1_charge_density susceptibility kernel
###########################


"""
Calculates the B1_charge_density susceptibility kernel 
"""
@inline Base.@propagate_inbounds function t0_diag_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10}, G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return t0_diag_B1_charge_density_kernel(mc, model, isks, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function t0_diag_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10},  G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return t0_diag_B1_charge_density_kernel(mc, model, isks, (G, G, G, G), flv, field)
end


@inline Base.@propagate_inbounds function t0_diag_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField  )
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 


    return nothing
end

@inline Base.@propagate_inbounds function t0_diag_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::Union{Abstract_DiscreteMBF1_X_symm, Cont_MBF1_X_symm, Discrete_MBF2_symm})
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 
    
    return (1/N)*real(2*((id - G0l.val[i, i])*Gl0.val[iMa1, iMa1] + (id - G0l.val[i, i])*Gl0.val[iMa2, iMa2] + (id - G0l.val[i, i])*Gl0.val[iPa1, iPa1] + 
    (id - G0l.val[i, i])*Gl0.val[iPa2, iPa2] + (id - G0l.val[i + N, i + N])*Gl0.val[iMa1 + N, iMa1 + N] + 
    (id - G0l.val[i + N, i + N])*Gl0.val[iMa2 + N, iMa2 + N] + (id - G0l.val[i + N, i + N])*Gl0.val[iPa1 + N, iPa1 + N] + 
    (id - G0l.val[i + N, i + N])*Gl0.val[iPa2 + N, iPa2 + N]))
end

###########################
### Zero-hopping, diagonal-only, Test 1, (t=0) B1_charge_density susceptibility kernel
###########################


"""
Calculates the B1_charge_density susceptibility kernel 
"""
@inline Base.@propagate_inbounds function t0_diag_test1_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10}, G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return t0_diag_test1_B1_charge_density_kernel(mc, model, isks, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function t0_diag_test1_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10},  G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return t0_diag_test1_B1_charge_density_kernel(mc, model, isks, (G, G, G, G), flv, field)
end


@inline Base.@propagate_inbounds function t0_diag_test1_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField  )
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 


    return nothing
end

@inline Base.@propagate_inbounds function t0_diag_test1_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::Union{Abstract_DiscreteMBF1_X_symm, Cont_MBF1_X_symm, Discrete_MBF2_symm})
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return (1/N)*real(-8*G0l.val[i, i]*Gl0.val[iPa1, iPa1] - 8*G0l.val[i + N, i + N]*Gl0.val[iPa1 + N, iPa1 + N])
end







###########################
### Local displaced charge_density susceptibility 
###########################
function local_displaced_charge_density_measurement(
    dqmc::DQMC, model::Model,  greens_iterator; 
        lattice_iterator = OneSite_Displaced(1, 2,-1), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = local_displaced_charge_density_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end


local_displaced_charge_density_correlation(args...; kwargs...) = local_displaced_charge_density_measurement(args..., Greens(); kwargs...)

local_displaced_charge_density_susceptibility(mc, args...; kwargs...) = local_displaced_charge_density_measurement(mc, args..., TimeIntegral(mc); kwargs...)


###########################
### local_displaced_charge_density susceptibility kernel
###########################
@inline Base.@propagate_inbounds function local_displaced_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    iiPt::NTuple{2}, G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return local_displaced_charge_density_kernel(mc, model, iiPt, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function local_displaced_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    iiPt::NTuple{2},  G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return local_displaced_charge_density_kernel(mc, model, iiPt, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function local_displaced_charge_density_kernel(
        mc, model::TwoBandModel, iiPt::NTuple{2}, packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField  )
    i, iPt = iiPt   # i, i+t
    μ, μp =flv
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    return (id*I[μ, μp]-G0l.val[i+(μ-1)*N, i+(μp-1)*N]) * 
        Gl0.val[iPt+(μp-1)*N, iPt+(μ-1)*N]
end


###########################
### lattice iterator for local_displaced charge density susceptibility
###########################

struct OneSite_Displaced <: DirectLatticeIterator 
    site::Int
    dir::Int
    l0::Int
end
OneSite_Displaced(::MonteCarloFlavor) = OneSite_Displaced(1, 2, -1)
output_size(::OneSite_Displaced, l::Lattice) = (1, )

_length(::OneSite_Displaced, l::Lattice) = 1

function apply!(
    temp::Array, iter::OneSite_Displaced, measurement, mc::DQMC, 
    packed_greens, weight = 1.0
    )
    G00, G0l, Gl0, Gll = packed_greens
    curr_l=G0l.l
    i=iter.site
    t=iter.dir
    l0=iter.l0

    @timeit_debug "apply!(::EachSitePair_B1, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        srcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        iPt=srcdir2trg[i, t]           # i + t
        if l0 ==-1 
            temp[1] += weight *measurement.kernel(mc, mc.model, 
            (i, iPt), packed_greens, (1, 1))
        elseif l0 ==curr_l
            temp[1] += measurement.kernel(mc, mc.model, 
            (i, iPt), packed_greens, (1, 1))
        end
                
    end
    return 
end

@inline function finalize_temp!(::OneSite_Displaced, m, mc)
    nothing
end
@inline function commit!(::OneSite_Displaced, m) 
    push!(m.observable, m.temp[1])
end


