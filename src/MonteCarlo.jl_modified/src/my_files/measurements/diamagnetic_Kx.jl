
###########################
### Kₓ onsite kinetic energy measurement, 
### a.k.a. diamagnetic contribution to phase stiffness
###########################
function kx_dia_measurement(
    dqmc::DQMC, model::Model; 
        greens_iterator=Greens(),
        lattice_iterator = EachWeightedBond(dqmc), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = kx_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end





@inline Base.@propagate_inbounds function kx_kernel(mc::DQMC, model::TwoBandModel, st::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return kx_kernel(mc, model, st, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function kx_kernel(mc::DQMC, model::TwoBandModel, st::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return kx_kernel(mc, model, st, (G, G, G, G), flv, field)
end
@inline Base.@propagate_inbounds function kx_kernel(mc, model, st::NTuple{2}, packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField)
    r, rPb =st
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(mc))
    T = mc.stack.hopping_matrix

    return G00.val[rPb, r]*T[r, rPb] + G00.val[N + rPb, N + r]*T[N + r, N + rPb] + 
        G00.val[2N + rPb, 2N + r]*T[2N + r, 2N + rPb] + G00.val[3N + rPb, 3N + r]*T[3N + r, 3N + rPb] + 
        G00.val[r, rPb]*T[rPb, r] + G00.val[N + r, N + rPb]*T[N + rPb, N + r] + 
        G00.val[2N + r, 2N + rPb]*T[2N + rPb, 2N + r] + G00.val[3N + r, 3N + rPb]*T[3N + rPb, 3N + r]
end
@inline Base.@propagate_inbounds function kx_kernel(mc, model, st::NTuple{2}, packed_greens::_GM4{<: Matrix}, 
    flv, ::Union{Discrete_MBF1_X_symm, Discrete_MBF2_symm})
    r, rPb =st
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(mc))
    T = mc.stack.hopping_matrix

    return 2*real(G00.val[rPb, r]*T[r, rPb] + G00.val[N + rPb, N + r]*T[N + r, N + rPb] + 
        G00.val[r, rPb]*T[rPb, r] + G00.val[N + r, N + rPb]*T[N + rPb, N + r])
end
