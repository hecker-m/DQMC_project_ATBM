

###########################
### nematic susceptibility 
###########################
function nematic_measurement(
    dqmc::DQMC, model::Model,  greens_iterator; 
        lattice_iterator = EachDoubleSitePairByDistance(), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = full_nem_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end

"""
    nematic_correlation(mc, model; kwargs...)

Generates an equal-time nematic correlation measurement. Note that the result needs to be added to the simulation 
via `mc[:name] = result`.

Note that internally in all three cases [Ising, XY, Heisenberg],
the expression is summed over all magnetization directions ∑_{ζ,ζ` ∈[x,y,z]}…

## Optional Keyword Arguments
- kwargs from `DQMCMeasurement`
"""
nematic_correlation(args...; kwargs...) = nematic_measurement(args..., Greens(); kwargs...)

"""
    nematic_susceptibility(mc, model; kwargs...)

Generates an time-integrated nematic susceptibility measurement. Note that the result needs to be added to the 
simulation via `mc[:name] = result`.

Note that internally in all three cases [Ising, XY, Heisenberg],
the expression is summed over all magnetization directions ∑_{ζ,ζ` ∈[x,y,z]}…

## Optional Keyword Arguments
- kwargs from `DQMCMeasurement`
"""
nematic_susceptibility(mc, args...; kwargs...) = nematic_measurement(mc, args..., TimeIntegral(mc); kwargs...)


###########################
### nematic susceptibility kernel
###########################


"""
Calculates the nematic susceptibility kernel 
"""
@inline Base.@propagate_inbounds function full_nem_kernel(mc::DQMC, model::TwoBandModel, klkPlP::NTuple{4}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_nem_kernel(mc, model, klkPlP, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_nem_kernel(mc::DQMC, model::TwoBandModel, klkPlP::NTuple{4}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_nem_kernel(mc, model, klkPlP, (G, G, G, G), flv, field)
end



@inline Base.@propagate_inbounds function full_nem_kernel(
        mc, model::TwoBandModel, klkPlP::NTuple{4}, packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField  )
    k, l, kPi, lPj = klkPlP   #k, l, k+i, l+j
	  G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    tcoscos1=(-(Gll.val[kPi, k + 3N]*Gll.val[k + N, kPi + 2N]) - Gll.val[k, kPi + 3N]*Gll.val[kPi + N, k + 2N] + Gll.val[k + N, kPi + 3N]*(2*Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
    Gll.val[k, kPi + 2N]*(Gll.val[kPi, k + 2N] + 2*Gll.val[kPi + N, k + 3N]) - Gll.val[kPi + N, k]*Gll.val[k + 2N, kPi + 3N] + 
    (Gll.val[k, kPi] + 2*Gll.val[k + N, kPi + N])*Gll.val[kPi + 2N, k + 2N] - Gll.val[k + N, kPi]*Gll.val[kPi + 2N, k + 3N] - 
    2*(Gll.val[kPi, kPi + 3N] + Gll.val[kPi + 2N, kPi + N])*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) - Gll.val[kPi + 2N, k + N]*Gll.val[k + 3N, kPi] + 
    Gll.val[kPi + 2N, k]*(Gll.val[k + 2N, kPi] + 2*Gll.val[k + 3N, kPi + N]) - Gll.val[kPi, k + N]*Gll.val[k + 3N, kPi + 2N] - Gll.val[k + 2N, kPi + N]*Gll.val[kPi + 3N, k] - 
    2*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N])*(Gll.val[kPi + N, kPi + 2N] + Gll.val[kPi + 3N, kPi]) + (2*Gll.val[k + 2N, kPi] + Gll.val[k + 3N, kPi + N])*Gll.val[kPi + 3N, k + N] - 
    (Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N])*(Gll.val[kPi, kPi + 2N] - Gll.val[kPi + N, kPi + 3N] + Gll.val[kPi + 2N, kPi] - 
      Gll.val[kPi + 3N, kPi + N]) - Gll.val[k, kPi + N]*Gll.val[kPi + 3N, k + 2N] + (2*Gll.val[k, kPi] + Gll.val[k + N, kPi + N])*Gll.val[kPi + 3N, k + 3N] + 
    Gll.val[kPi + N, k + N]*(2*Gll.val[k + 2N, kPi + 2N] + Gll.val[k + 3N, kPi + 3N] - 3*I[k, kPi]) + 
    Gll.val[kPi, k]*(Gll.val[k + 2N, kPi + 2N] + 2*Gll.val[k + 3N, kPi + 3N] - 3*I[k, kPi]) - 3*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N])*I[k, kPi])*
   (-(G00.val[lPj, l + 3N]*G00.val[l + N, lPj + 2N]) - G00.val[l, lPj + 3N]*G00.val[lPj + N, l + 2N] + G00.val[l + N, lPj + 3N]*(2*G00.val[lPj, l + 2N] + G00.val[lPj + N, l + 3N]) + 
    G00.val[l, lPj + 2N]*(G00.val[lPj, l + 2N] + 2*G00.val[lPj + N, l + 3N]) - G00.val[lPj + N, l]*G00.val[l + 2N, lPj + 3N] + 
    (G00.val[l, lPj] + 2*G00.val[l + N, lPj + N])*G00.val[lPj + 2N, l + 2N] - G00.val[l + N, lPj]*G00.val[lPj + 2N, l + 3N] - 
    2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*(G00.val[l + N, l + 2N] + G00.val[l + 3N, l]) - G00.val[lPj + 2N, l + N]*G00.val[l + 3N, lPj] + 
    G00.val[lPj + 2N, l]*(G00.val[l + 2N, lPj] + 2*G00.val[l + 3N, lPj + N]) - G00.val[lPj, l + N]*G00.val[l + 3N, lPj + 2N] - G00.val[l + 2N, lPj + N]*G00.val[lPj + 3N, l] - 
    2*(G00.val[l, l + 3N] + G00.val[l + 2N, l + N])*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj]) + (2*G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G00.val[lPj + 3N, l + N] - 
    (G00.val[l, l + 2N] - G00.val[l + N, l + 3N] + G00.val[l + 2N, l] - G00.val[l + 3N, l + N])*(G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - 
      G00.val[lPj + 3N, lPj + N]) - G00.val[l, lPj + N]*G00.val[lPj + 3N, l + 2N] + (2*G00.val[l, lPj] + G00.val[l + N, lPj + N])*G00.val[lPj + 3N, l + 3N] + 
    G00.val[lPj + N, l + N]*(2*G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N] - 3*I[l, lPj]) + 
    G00.val[lPj, l]*(G00.val[l + 2N, lPj + 2N] + 2*G00.val[l + 3N, lPj + 3N] - 3*I[l, lPj]) - 3*(G00.val[lPj + 2N, l + 2N] + G00.val[lPj + 3N, l + 3N])*I[l, lPj]);


    tcoscos2=(G0l.val[l + 3N, kPi + 2N]*G0l.val[lPj + 3N, k + 2N] - G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + 2N])*
    (Gl0.val[k, lPj + N]*Gl0.val[kPi, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi, lPj + N]) + (G0l.val[l, kPi + 2N]*G0l.val[lPj, k + 2N] - G0l.val[l, k + 2N]*G0l.val[lPj, kPi + 2N])*
    (Gl0.val[k, lPj + 2N]*Gl0.val[kPi, l + 2N] - Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + 2N]) + 
   2*(G0l.val[l, kPi + 2N]*G0l.val[lPj + 3N, k + 2N] - G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + 2N])*(-2*Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj] + Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + N] - 
     Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 2N] + 2*Gl0.val[k, lPj]*Gl0.val[kPi, l + 3N]) + 
   2*(-(G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 3N, k + 2N]) + G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + 2N])*
    (Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 3N]) - 
   2*(G0l.val[lPj, kPi + 2N]*G0l.val[l + N, k + 2N] - G0l.val[lPj, k + 2N]*G0l.val[l + N, kPi + 2N] + 2*G0l.val[l, kPi + 2N]*G0l.val[lPj + N, k + 2N] - 
     2*G0l.val[l, k + 2N]*G0l.val[lPj + N, kPi + 2N])*(Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi, l + 3N]) + 
   (G0l.val[l + N, kPi + 2N]*G0l.val[lPj + N, k + 2N] - G0l.val[l + N, k + 2N]*G0l.val[lPj + N, kPi + 2N])*(Gl0.val[k, lPj + 3N]*Gl0.val[kPi, l + 3N] - 
     Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + 3N]) + (G0l.val[l + 2N, kPi + 3N]*G0l.val[lPj + 2N, k + 3N] - G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi + 3N])*
    (Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l] - Gl0.val[k + N, l]*Gl0.val[kPi + N, lPj]) + 2*(-(G0l.val[l, kPi + 3N]*G0l.val[lPj + 2N, k + 3N]) + G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi + 3N])*
    (Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 2N]) - 2*(G0l.val[l, k + 3N]*G0l.val[lPj, kPi + 2N] - G0l.val[l, kPi + 2N]*G0l.val[lPj, k + 3N])*
    (Gl0.val[kPi, lPj + 2N]*Gl0.val[k + N, l + 2N] - Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + 2N] + 2*Gl0.val[k, lPj + 2N]*Gl0.val[kPi + N, l + 2N] - 
     2*Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj + 2N]) + (G0l.val[l, kPi + 3N]*G0l.val[lPj, k + 3N] - G0l.val[l, k + 3N]*G0l.val[lPj, kPi + 3N])*
    (Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + N, l + 2N] - Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + 2N]) + 
   4*(G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi + 2N] - G0l.val[lPj, kPi + 2N]*G0l.val[l + N, k + 3N] + 2*G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi + 2N] - 
     2*G0l.val[l, kPi + 2N]*G0l.val[lPj + N, k + 3N])*(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + 2N] - Gl0.val[kPi, lPj + 2N]*Gl0.val[k + N, l + 3N] + 
     2*Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj + 2N] - 2*Gl0.val[k, lPj + 2N]*Gl0.val[kPi + N, l + 3N]) + 
   2*(G0l.val[l + N, kPi + 3N]*G0l.val[lPj + 2N, k + 3N] - G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + 3N])*
    (Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj] - 2*Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + N] + 2*Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 2N] - 
     Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 3N]) - 2*(G0l.val[lPj, kPi + 3N]*G0l.val[l + N, k + 3N] - G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi + 3N] + 
     2*G0l.val[l, kPi + 3N]*G0l.val[lPj + N, k + 3N] - 2*G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi + 3N])*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + 2N] - 
     Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + N, l + 3N]) - 2*(G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi + 2N] - G0l.val[l + N, kPi + 2N]*G0l.val[lPj + N, k + 3N])*
    (Gl0.val[kPi, lPj + 3N]*Gl0.val[k + N, l + 3N] - Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + 3N] + 2*Gl0.val[k, lPj + 3N]*Gl0.val[kPi + N, l + 3N] - 
     2*Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj + 3N]) + (G0l.val[l + N, kPi + 3N]*G0l.val[lPj + N, k + 3N] - G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi + 3N])*
    (Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + N, l + 3N] - Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + 3N]) + 
   (G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k] - G0l.val[l + 2N, k]*G0l.val[lPj + 2N, kPi])*(Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l] - Gl0.val[k + 2N, l]*Gl0.val[kPi + 2N, lPj]) - 
   2*(G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k] - G0l.val[lPj + 2N, k]*G0l.val[l + 3N, kPi] + 2*G0l.val[l + 2N, kPi]*G0l.val[lPj + 3N, k] - 2*G0l.val[l + 2N, k]*G0l.val[lPj + 3N, kPi])*
    (Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + N]) + 
   2*(G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k + 2N])*(-(Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + N]) + 
     Gl0.val[k, l + N]*Gl0.val[kPi + 2N, lPj + N]) + (G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k] - G0l.val[l + 3N, k]*G0l.val[lPj + 3N, kPi])*
    (Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + N] - Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 2N, lPj + N]) + 
   4*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k + 2N])*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - 
     Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 3N]) + 2*(-(G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k]) + G0l.val[l + N, k]*G0l.val[lPj + 3N, kPi])*
    (Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + 3N]) + 
   2*(G0l.val[l + N, k + 2N]*G0l.val[lPj + N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 2N])*(-(Gl0.val[k, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N]) + 
     Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N]) + (G0l.val[l + N, kPi]*G0l.val[lPj + N, k] - G0l.val[l + N, k]*G0l.val[lPj + N, kPi])*
    (Gl0.val[k + 2N, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N]) + 
   2*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + 3N])*(Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l] - Gl0.val[k + N, l]*Gl0.val[kPi + 2N, lPj] - 
     2*Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l] + 2*Gl0.val[k, l]*Gl0.val[kPi + 3N, lPj]) + 
   2*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 2N, k + 3N])*
    (-(Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l]) + Gl0.val[k + N, l]*Gl0.val[kPi + 3N, lPj]) - 2*(G0l.val[l + 2N, k + N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + N])*
    (Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l] - Gl0.val[kPi + 2N, l]*Gl0.val[k + 3N, lPj] + 2*Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l] - 2*Gl0.val[k + 2N, l]*Gl0.val[kPi + 3N, lPj]) + 
   (G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 2N, k + N] - G0l.val[l + 2N, k + N]*G0l.val[lPj + 2N, kPi + N])*(Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l] - 
     Gl0.val[k + 3N, l]*Gl0.val[kPi + 3N, lPj]) + 4*(G0l.val[lPj + 2N, k + N]*G0l.val[l + 3N, kPi] - G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k + N] + 
     2*G0l.val[l + 2N, k + N]*G0l.val[lPj + 3N, kPi] - 2*G0l.val[l + 2N, kPi]*G0l.val[lPj + 3N, k + N])*(Gl0.val[kPi + 2N, l + N]*Gl0.val[k + 3N, lPj] - 
     Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + N] + 2*Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 3N, lPj] - 2*Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + N]) - 
   2*(G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + N] - G0l.val[lPj + 2N, k + N]*G0l.val[l + 3N, kPi + N] + 2*G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 3N, k + N] - 
     2*G0l.val[l + 2N, k + N]*G0l.val[lPj + 3N, kPi + N])*(Gl0.val[k + 3N, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + N]) - 
   2*(G0l.val[l + 3N, k + N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k + N])*(Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + N] - 
     Gl0.val[kPi + 2N, l + N]*Gl0.val[k + 3N, lPj + N] + 2*Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + N] - 2*Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 3N, lPj + N]) + 
   (G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + N] - G0l.val[l + 3N, k + N]*G0l.val[lPj + 3N, kPi + N])*(Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + N] - 
     Gl0.val[k + 3N, l + N]*Gl0.val[kPi + 3N, lPj + N]) + 4*(G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj + 2N, k + 3N])*
    (Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + 2N]) + 
   2*(-(G0l.val[l, kPi + N]*G0l.val[lPj + 2N, k + N]) + G0l.val[l, k + N]*G0l.val[lPj + 2N, kPi + N])*(Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj] - 
     Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + 2N]) + 2*(G0l.val[l, k + 3N]*G0l.val[lPj, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj, k + 3N])*
    (-(Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N]) + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N]) + 
   (G0l.val[l, kPi + N]*G0l.val[lPj, k + N] - G0l.val[l, k + N]*G0l.val[lPj, kPi + N])*(Gl0.val[k + 3N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N] - 
     Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N]) + 4*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + 3N])*
    (Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] - 2*Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] + 2*Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N] - 
     Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N] - 2*Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] + 4*Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - 
     4*Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N] + 2*Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N]) - 
   4*(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + N] - Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 3N] + 2*Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - 
     2*Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(-(G0l.val[l + N, k + N]*G0l.val[lPj + 3N, kPi]) + G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k + N] + 
     id*G0l.val[lPj + 3N, kPi]*I[l, k]) + 2*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + 2N])*
    (-2*G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k] + 2*G0l.val[l + N, k]*G0l.val[lPj + 2N, kPi] + G0l.val[lPj + 3N, kPi]*(-G0l.val[l, k] + id*I[l, k]) + 
     G0l.val[lPj + 3N, k]*(G0l.val[l, kPi] - id*I[l, kPi])) - 2*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + 3N])*
    (-(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k]) + G0l.val[l + N, k]*G0l.val[lPj + 2N, kPi] - 2*G0l.val[lPj + 3N, kPi]*(G0l.val[l, k] - id*I[l, k]) + 
     2*G0l.val[lPj + 3N, k]*(G0l.val[l, kPi] - id*I[l, kPi])) + 2*(Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*
    (G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + N] - G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi + N] + 2*G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) - 
     2*G0l.val[lPj + 2N, k + N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi])) + 2*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + 2N])*
    (G0l.val[lPj + 2N, kPi]*(G0l.val[l, k] - id*I[l, k]) + G0l.val[lPj + 2N, k]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 
   4*(-(Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj]) + Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 2N] - 2*Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj] + 
     2*Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + 2N])*(G0l.val[l, k + N]*G0l.val[lPj + 2N, kPi] + G0l.val[lPj + 2N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 
   4*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 2N] - 2*Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj] + 2*Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 2N])*
    (G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi] + G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
   4*(2*Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + N] + Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 2N] - 2*Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 3N])*
    (G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi] + G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 
   2*(Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + 3N])*(2*G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + N] - 
     2*G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + 
     G0l.val[lPj + 2N, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 2*(Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*
    (G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[lPj + 3N, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 
   4*((Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj] - 2*Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + N] + 2*Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 2N] - 
       Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 3N] + 2*Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj] - 4*Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] + 
       4*Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + 2N] - 2*Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + 3N])*(-(G0l.val[l + N, k + N]*G0l.val[lPj + 2N, kPi]) + 
       G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + N] + id*G0l.val[lPj + 2N, kPi]*I[l, k]) + 
     (2*Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj] - Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + N] + Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 2N] - 
       2*Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 3N] + 4*Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj] - 2*Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] + 
       2*Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + 2N] - 4*Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi] + 
       G0l.val[lPj + 3N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi]))) + 4*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*
    (2*G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + N] - 2*G0l.val[lPj + 2N, k + 3N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) + 
     G0l.val[l, kPi + N]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k])) + 4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*
    (2*G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[lPj + 3N, k + 2N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) - 
     2*G0l.val[l + N, kPi]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k])) + 2*(-(Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + N]) + Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj + N])*
    (-(G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + N]) + G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + 2N] + 2*G0l.val[lPj + 3N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) - 
     2*G0l.val[l + 3N, kPi]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k])) - 
   2*(Gl0.val[kPi + 2N, lPj + 3N]*Gl0.val[k + 3N, l + 3N] - Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + 3N] + 2*Gl0.val[k + 2N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N] - 
     2*Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*(G0l.val[lPj + N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[l + N, kPi]*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) + 
   4*(2*Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj] - 2*Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 2N] - Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 2N])*
    (G0l.val[l, k + 2N]*G0l.val[lPj + 2N, kPi + N] + G0l.val[l, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 
   4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj] - 2*Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + N] + 2*Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 2N] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 3N])*
    (G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi] + G0l.val[l + N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 
   2*(-(Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l]) + Gl0.val[k, l]*Gl0.val[kPi + 2N, lPj])*(G0l.val[lPj + 2N, kPi]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
     G0l.val[l + 2N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 2*(-2*Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l] + 2*Gl0.val[k + N, l]*Gl0.val[kPi + 2N, lPj] + 
     Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l] - Gl0.val[k, l]*Gl0.val[kPi + 3N, lPj])*(G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
     G0l.val[l + 2N, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 4*(Gl0.val[k, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + N])*
    (G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k + 2N] + 2*G0l.val[l + 2N, kPi]*G0l.val[lPj + 3N, k + 2N] - 2*G0l.val[lPj + 3N, kPi]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
     G0l.val[l + 3N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 4*(Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + N])*
    (2*G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi] + G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + 2N] + 4*G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi] + 
     2*G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 3N, k + 2N] - 2*G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) - 
     2*G0l.val[lPj + 2N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) - 4*G0l.val[l + 2N, kPi]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k]) + 
     G0l.val[l + 3N, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 4*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 2N])*
    (G0l.val[l, k + 2N]*G0l.val[lPj + 2N, kPi] + (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 
   4*(2*Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj] - Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + N] + Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 2N] - 2*Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 3N] + 
     4*Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj] - 2*Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj + N] + 2*Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + 2N] - 4*Gl0.val[k, lPj]*Gl0.val[kPi + N, l + 3N])*
    (G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
   4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi] - 
     2*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 
     G0l.val[l + N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 4*(-(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + N]) + Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 3N] - 
     2*Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj + N] + 2*Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + 
     G0l.val[l + N, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 2*(-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj + N])*
    (-2*G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + N] + 2*G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + 2N] + G0l.val[lPj + 3N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + 
     G0l.val[l + 3N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 2*(-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + 3N, lPj + N])*
    (G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + G0l.val[l + 3N, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
   2*(Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + N] - Gl0.val[kPi, l + N]*Gl0.val[k + N, lPj + N] + 2*Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + N] - 2*Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj + N])*
    (G0l.val[lPj + 3N, kPi + 2N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + G0l.val[l + 3N, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
   4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi + N] + 
     (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
   4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + N] + 
     G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) - 
     2*(G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[l, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]))) + 
   4*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + N])*(-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi]) + 
     G0l.val[lPj + 2N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) - 2*(G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + 2N] + G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi] + 
       2*G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 3N, k + 2N] - 2*G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
       G0l.val[l + 3N, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l + 2N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]))) - 
   4*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + N])*(-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi + N]) + 
     G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) - 2*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi + N] + 
       G0l.val[l + 2N, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]))) + 
   4*((G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + 2N])*(4*Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] - 
       2*Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] + 2*Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N] - 4*Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N] - 
       2*Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N] + 
       2*Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N]) + (-2*Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] + 4*Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - 
       4*Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N] + 2*Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N] + Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] - 
       2*Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] + 2*Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N])*
      (G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi + N] + (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 
     (-2*Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N] + 
       2*Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N] + 4*Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] - 2*Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] + 
       2*Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N] - 4*Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi] + 
       (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]))) + 
   2*(-(Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N])*(G0l.val[lPj, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
     G0l.val[l, k + 2N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 2*(-(Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N])*
    (-2*G0l.val[l, k + 2N]*G0l.val[lPj, kPi + N] + 2*G0l.val[l, kPi + N]*G0l.val[lPj, k + 2N] + G0l.val[lPj, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
     G0l.val[l, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 2*(-(Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N]) + Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N])*
    (-(G0l.val[l, k + 2N]*G0l.val[lPj, kPi + N]) + G0l.val[l, kPi + N]*G0l.val[lPj, k + 2N] - 2*G0l.val[lPj, k + 3N]*(G0l.val[l, kPi] - id*I[l, kPi]) + 
     2*G0l.val[l, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 2*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
    (-2*G0l.val[lPj + N, kPi]*(G0l.val[l, k] - id*I[l, k]) + 2*G0l.val[lPj + N, k]*(G0l.val[l, kPi] - id*I[l, kPi]) + G0l.val[l + N, kPi]*(-G0l.val[lPj, k] + id*I[lPj, k]) + 
     G0l.val[l + N, k]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
    (-(G0l.val[lPj, k + 2N]*G0l.val[l + N, kPi]) - 2*(G0l.val[l, k + 2N]*G0l.val[lPj + N, kPi] + G0l.val[lPj + N, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
     G0l.val[l + N, k + 2N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 2*(Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + 3N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
    (G0l.val[lPj, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[lPj, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + 
     2*G0l.val[l, kPi + N]*(G0l.val[lPj + N, k + N] - id*I[lPj, k]) - 2*G0l.val[l, k + N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
   4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*(2*G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi] + 
     G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 2N] + 4*G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi] + 2*G0l.val[l, kPi + N]*G0l.val[lPj + N, k + 2N] - 
     4*G0l.val[lPj + N, k + 3N]*(G0l.val[l, kPi] - id*I[l, kPi]) + G0l.val[lPj, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) - 
     2*G0l.val[l + N, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi]) - 2*G0l.val[l, k + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) - 
   4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
    (G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 3N] + 2*G0l.val[l, kPi + N]*G0l.val[lPj + N, k + 3N] + G0l.val[lPj, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) - 
     2*G0l.val[l, k + 3N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) - 2*(-(Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N]) + Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N])*
    (G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 3N] + 2*G0l.val[lPj + N, k + 2N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) - 
     2*G0l.val[l + N, k + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
   2*(-(Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*(G0l.val[lPj + N, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + 
     G0l.val[l + N, k + 3N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 2*(Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + 2N])*
    (G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
   2*(Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj] - 2*Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + N] + 2*Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 2N] - Gl0.val[k, lPj]*Gl0.val[kPi, l + 3N])*
    (G0l.val[l + N, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l + N, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
   2*(Gl0.val[kPi, lPj]*Gl0.val[k + N, l] - Gl0.val[kPi, l]*Gl0.val[k + N, lPj] + 2*Gl0.val[k, lPj]*Gl0.val[kPi + N, l] - 2*Gl0.val[k, l]*Gl0.val[kPi + N, lPj])*
    (G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + 2N, kPi + 2N] + id*I[l, kPi]) + G0l.val[l + 2N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
   2*(Gl0.val[k, l + N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + N])*(-2*G0l.val[lPj + 3N, kPi + 2N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
     2*G0l.val[lPj + 3N, k + 2N]*(G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi]) + G0l.val[l + 3N, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
     G0l.val[l + 3N, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
   2*(2*Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + N] + Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 2N] - 
     2*Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 3N])*(G0l.val[l, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + 
     G0l.val[l, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 2*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 3N])*
    (G0l.val[l + N, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) - 
   2*(Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + N])*(G0l.val[lPj + 2N, kPi + 3N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + 
     G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + 3N, kPi + 3N] + id*I[l, kPi]) + 2*G0l.val[l + 2N, kPi + 3N]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k]) - 
     2*G0l.val[l + 2N, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 
   2*(Gl0.val[kPi + 2N, lPj + 2N]*Gl0.val[k + 3N, l + 2N] - Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + 2N] + 2*Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N] - 
     2*Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N])*(G0l.val[lPj, k + N]*(G0l.val[l, kPi] - id*I[l, kPi]) + G0l.val[l, k + N]*(-G0l.val[lPj, kPi] + id*I[lPj, kPi])) + 
   (-(Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N])*
    ((G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj, k] + id*I[lPj, k]) + (-G0l.val[l, k] + id*I[l, k])*(-G0l.val[lPj, kPi] + id*I[lPj, kPi])) + 
   4*(-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + 2N]) + Gl0.val[kPi + 2N, lPj + 2N]*Gl0.val[k + 3N, l + 3N] - 2*Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] + 
     2*Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*(-(G0l.val[lPj, k + N]*G0l.val[l + N, kPi]) - 2*G0l.val[l, k + N]*G0l.val[lPj + N, kPi] - 
     2*(G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k]) + (-G0l.val[l + N, k + N] + id*I[l, k])*(-G0l.val[lPj, kPi] + id*I[lPj, kPi])) + 
   2*(-(Gl0.val[k, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*(2*G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi] - 
     2*G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 3N] + G0l.val[lPj + N, k + 2N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) + 
     G0l.val[l + N, k + 2N]*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi])) + 
   (-(Gl0.val[k + 3N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*
    ((G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k]) + (-G0l.val[l + N, k + N] + id*I[l, k])*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi])) - 
   4*(Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj] - Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 2N] + 2*Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj] - 2*Gl0.val[k, lPj]*Gl0.val[kPi + N, l + 2N])*
    (G0l.val[l, kPi + 2N]*G0l.val[lPj + 2N, k + 3N] + G0l.val[l, k + 3N]*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi])) + 
   4*(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj] - 2*Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + N] + 2*Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 2N] - Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 3N] + 
     2*Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj] - 4*Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj + N] + 4*Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + 2N] - 2*Gl0.val[k, lPj]*Gl0.val[kPi + N, l + 3N])*
    (G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 2N, k + 3N] + G0l.val[l + N, k + 3N]*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi])) + 
   (-(Gl0.val[k, lPj]*Gl0.val[kPi, l]) + Gl0.val[k, l]*Gl0.val[kPi, lPj])*((G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
     (-G0l.val[l + 2N, k + 2N] + id*I[l, k])*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi])) + 
   4*(-(Gl0.val[kPi, l + N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + N] - 2*Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj] + 2*Gl0.val[k, lPj]*Gl0.val[kPi + N, l + N])*
    (-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi + 2N]) - 2*G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] - 
     2*(G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + (-G0l.val[l + 3N, k + 3N] + id*I[l, k])*
      (-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi])) + (-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj + N])*
    ((G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + (-G0l.val[l + 3N, k + 3N] + id*I[l, k])*
      (-G0l.val[lPj + 3N, kPi + 3N] + id*I[lPj, kPi])) + 4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
    (-(G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi]) - 
     2*(G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 2N] + G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi] + 2*G0l.val[l, kPi + N]*G0l.val[lPj + N, k + 2N] + 
       G0l.val[lPj + N, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + G0l.val[lPj, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) - 
       2*G0l.val[l, k + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])));

    tcoscos3=(-2*Gl0.val[kPi + 2N, l + 2N]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + Gl0.val[kPi + 3N, l + 2N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - 
    Gll.val[k + 3N, k + N]) - Gl0.val[k + 2N, l + 2N]*Gll.val[kPi + 3N, k] + Gl0.val[k + 3N, l + 2N]*(2*Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
  Gl0.val[k, l + 2N]*Gll.val[kPi + 3N, k + 2N] + Gl0.val[k + N, l + 2N]*(2*Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
 ((G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l, kPi + N] - 
  (G00.val[l, lPj + 2N] + 2*G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi + N] + G00.val[l, lPj + 3N]*G0l.val[lPj + N, kPi + N] - 
  (G00.val[l, lPj] + 2*G00.val[l + N, lPj + N])*G0l.val[lPj + 2N, kPi + N] + G00.val[l, lPj + N]*G0l.val[lPj + 3N, kPi + N] + 
  2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) - id*G00.val[l, lPj + 3N]*I[lPj, kPi]) + 
(-2*Gl0.val[kPi + 3N, l + 3N]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) - Gl0.val[k + 3N, l + 3N]*Gll.val[kPi + 2N, k + N] - Gl0.val[k + N, l + 3N]*Gll.val[kPi + 2N, k + 3N] + 
  Gl0.val[kPi + 2N, l + 3N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]) + 
  Gl0.val[k + 2N, l + 3N]*(Gll.val[kPi + 2N, k] + 2*Gll.val[kPi + 3N, k + N]) + Gl0.val[k, l + 3N]*(Gll.val[kPi + 2N, k + 2N] + 2*Gll.val[kPi + 3N, k + 3N]))*
 (G00.val[l + N, lPj + 2N]*G0l.val[lPj, kPi] + (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*G0l.val[l + N, kPi] - 
  (2*G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi] + G00.val[l + N, lPj]*G0l.val[lPj + 2N, kPi] - (2*G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 3N, kPi] + 
  2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*(G0l.val[l, kPi] - id*I[l, kPi]) - id*G00.val[l + N, lPj + 2N]*I[lPj, kPi]) + 
(-2*Gl0.val[kPi + 2N, l + 3N]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + Gl0.val[kPi + 3N, l + 3N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - 
    Gll.val[k + 3N, k + N]) - Gl0.val[k + 2N, l + 3N]*Gll.val[kPi + 3N, k] + Gl0.val[k + 3N, l + 3N]*(2*Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
  Gl0.val[k, l + 3N]*Gll.val[kPi + 3N, k + 2N] + Gl0.val[k + N, l + 3N]*(2*Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
 (2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*G0l.val[l, kPi + N] + G00.val[l + N, lPj + 2N]*G0l.val[lPj, kPi + N] - 
  (2*G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi + N] + G00.val[l + N, lPj]*G0l.val[lPj + 2N, kPi + N] - 
  (2*G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 3N, kPi + N] - (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*
   (G0l.val[l + N, kPi + N] - id*I[l, kPi]) + id*(2*G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*I[lPj, kPi]) + 
(-2*Gl0.val[kPi + 3N, l + 2N]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) - Gl0.val[k + 3N, l + 2N]*Gll.val[kPi + 2N, k + N] - Gl0.val[k + N, l + 2N]*Gll.val[kPi + 2N, k + 3N] + 
  Gl0.val[kPi + 2N, l + 2N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]) + 
  Gl0.val[k + 2N, l + 2N]*(Gll.val[kPi + 2N, k] + 2*Gll.val[kPi + 3N, k + N]) + Gl0.val[k, l + 2N]*(Gll.val[kPi + 2N, k + 2N] + 2*Gll.val[kPi + 3N, k + 3N]))*
 (-((G00.val[l, lPj + 2N] + 2*G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi]) + 2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*G0l.val[l + N, kPi] + 
  G00.val[l, lPj + 3N]*G0l.val[lPj + N, kPi] - (G00.val[l, lPj] + 2*G00.val[l + N, lPj + N])*G0l.val[lPj + 2N, kPi] + G00.val[l, lPj + N]*G0l.val[lPj + 3N, kPi] + 
  (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*(G0l.val[l, kPi] - id*I[l, kPi]) + 
  id*(G00.val[l, lPj + 2N] + 2*G00.val[l + N, lPj + 3N])*I[lPj, kPi]) + (-(Gl0.val[k + 2N, l]*Gll.val[kPi + N, k]) + Gl0.val[k + 3N, l]*(2*Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
  Gl0.val[k, l]*Gll.val[kPi + N, k + 2N] + Gl0.val[k + N, l]*(2*Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) - 2*Gl0.val[kPi, l]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + 
  Gl0.val[kPi + N, l]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
 (-((G00.val[l + 2N, lPj + 2N] + 2*G00.val[l + 3N, lPj + 3N])*G0l.val[lPj, kPi + 3N]) + G00.val[l + 2N, lPj + 3N]*G0l.val[lPj + N, kPi + 3N] + 
  (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l + 2N, kPi + 3N] - 
  (G00.val[l + 2N, lPj] + 2*G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi + 3N] + G00.val[l + 2N, lPj + N]*G0l.val[lPj + 3N, kPi + 3N] + 
  2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*(G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi]) - id*G00.val[l + 2N, lPj + N]*I[lPj, kPi]) + 
(-(Gl0.val[k + 3N, l + N]*Gll.val[kPi, k + N]) - Gl0.val[k + N, l + N]*Gll.val[kPi, k + 3N] + Gl0.val[k + 2N, l + N]*(Gll.val[kPi, k] + 2*Gll.val[kPi + N, k + N]) + 
  Gl0.val[k, l + N]*(Gll.val[kPi, k + 2N] + 2*Gll.val[kPi + N, k + 3N]) - 2*Gl0.val[kPi + N, l + N]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) + 
  Gl0.val[kPi, l + N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]))*
 (G00.val[l + 3N, lPj + 2N]*G0l.val[lPj, kPi + 2N] - (2*G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj + N, kPi + 2N] + 
  G00.val[l + 3N, lPj]*G0l.val[lPj + 2N, kPi + 2N] + (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*
   G0l.val[l + 3N, kPi + 2N] - (2*G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi + 2N] + 
  2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*(G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi]) - id*G00.val[l + 3N, lPj]*I[lPj, kPi]) + 
(-(Gl0.val[k + 2N, l + N]*Gll.val[kPi + N, k]) + Gl0.val[k + 3N, l + N]*(2*Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - Gl0.val[k, l + N]*Gll.val[kPi + N, k + 2N] + 
  Gl0.val[k + N, l + N]*(2*Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) - 2*Gl0.val[kPi, l + N]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + 
  Gl0.val[kPi + N, l + N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
 (G00.val[l + 3N, lPj + 2N]*G0l.val[lPj, kPi + 3N] - (2*G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj + N, kPi + 3N] + 
  2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*G0l.val[l + 2N, kPi + 3N] + G00.val[l + 3N, lPj]*G0l.val[lPj + 2N, kPi + 3N] - 
  (2*G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi + 3N] - (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - 
    G00.val[lPj + 3N, lPj + N])*(G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi]) + id*(2*G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*I[lPj, kPi]) + 
(-(Gl0.val[k + 3N, l]*Gll.val[kPi, k + N]) - Gl0.val[k + N, l]*Gll.val[kPi, k + 3N] + Gl0.val[k + 2N, l]*(Gll.val[kPi, k] + 2*Gll.val[kPi + N, k + N]) + 
  Gl0.val[k, l]*(Gll.val[kPi, k + 2N] + 2*Gll.val[kPi + N, k + 3N]) - 2*Gl0.val[kPi + N, l]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) + 
  Gl0.val[kPi, l]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]))*
 (-((G00.val[l + 2N, lPj + 2N] + 2*G00.val[l + 3N, lPj + 3N])*G0l.val[lPj, kPi + 2N]) + G00.val[l + 2N, lPj + 3N]*G0l.val[lPj + N, kPi + 2N] - 
  (G00.val[l + 2N, lPj] + 2*G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi + 2N] + 2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*G0l.val[l + 3N, kPi + 2N] + 
  G00.val[l + 2N, lPj + N]*G0l.val[lPj + 3N, kPi + 2N] + (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*
   (G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi]) + id*(G00.val[l + 2N, lPj] + 2*G00.val[l + 3N, lPj + N])*I[lPj, kPi]) + 
(-2*Gl0.val[kPi + 3N, l]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) - Gl0.val[k + 3N, l]*Gll.val[kPi + 2N, k + N] - Gl0.val[k + N, l]*Gll.val[kPi + 2N, k + 3N] + 
  Gl0.val[kPi + 2N, l]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]) + 
  Gl0.val[k + 2N, l]*(Gll.val[kPi + 2N, k] + 2*Gll.val[kPi + 3N, k + N]) + Gl0.val[k, l]*(Gll.val[kPi + 2N, k + 2N] + 2*Gll.val[kPi + 3N, k + 3N]))*
 (G00.val[l + 2N, lPj + 3N]*G0l.val[lPj + N, kPi] + (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l + 2N, kPi] - 
  (G00.val[l + 2N, lPj] + 2*G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi] + 2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*G0l.val[l + 3N, kPi] + 
  G00.val[l + 2N, lPj + N]*G0l.val[lPj + 3N, kPi] - (G00.val[l + 2N, lPj + 2N] + 2*G00.val[l + 3N, lPj + 3N])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
(-2*Gl0.val[kPi + 2N, l + N]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + Gl0.val[kPi + 3N, l + N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - 
    Gll.val[k + 3N, k + N]) - Gl0.val[k + 2N, l + N]*Gll.val[kPi + 3N, k] + Gl0.val[k + 3N, l + N]*(2*Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
  Gl0.val[k, l + N]*Gll.val[kPi + 3N, k + 2N] + Gl0.val[k + N, l + N]*(2*Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
 (G00.val[l + 3N, lPj + 2N]*G0l.val[lPj, kPi + N] + 2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*G0l.val[l + 2N, kPi + N] + G00.val[l + 3N, lPj]*G0l.val[lPj + 2N, kPi + N] + 
  (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*G0l.val[l + 3N, kPi + N] - 
  (2*G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi + N] - (2*G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*
   (G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) - (Gl0.val[k + 3N, l + 3N]*Gll.val[kPi, k + N] + Gl0.val[k + N, l + 3N]*Gll.val[kPi, k + 3N] - 
  Gl0.val[k + 2N, l + 3N]*(Gll.val[kPi, k] + 2*Gll.val[kPi + N, k + N]) - Gl0.val[k, l + 3N]*(Gll.val[kPi, k + 2N] + 2*Gll.val[kPi + N, k + 3N]) + 
  2*Gl0.val[kPi + N, l + 3N]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) + Gl0.val[kPi, l + 3N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - 
    Gll.val[k + 3N, k + N]))*(2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*G0l.val[l, kPi + 2N] + G00.val[l + N, lPj + 2N]*G0l.val[lPj, kPi + 2N] + 
  (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*G0l.val[l + N, kPi + 2N] - 
  (2*G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi + 2N] - (2*G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 3N, kPi + 2N] + 
  G00.val[l + N, lPj]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + (-(Gl0.val[k + 3N, l + 2N]*Gll.val[kPi, k + N]) - Gl0.val[k + N, l + 2N]*Gll.val[kPi, k + 3N] + 
  Gl0.val[k + 2N, l + 2N]*(Gll.val[kPi, k] + 2*Gll.val[kPi + N, k + N]) + Gl0.val[k, l + 2N]*(Gll.val[kPi, k + 2N] + 2*Gll.val[kPi + N, k + 3N]) - 
  2*Gl0.val[kPi + N, l + 2N]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) + Gl0.val[kPi, l + 2N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + 
    Gll.val[k + 3N, k + N]))*((G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l, kPi + 2N] - 
  (G00.val[l, lPj + 2N] + 2*G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi + 2N] + 2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*G0l.val[l + N, kPi + 2N] + 
  G00.val[l, lPj + 3N]*G0l.val[lPj + N, kPi + 2N] + G00.val[l, lPj + N]*G0l.val[lPj + 3N, kPi + 2N] - (G00.val[l, lPj] + 2*G00.val[l + N, lPj + N])*
   (G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + (-(Gl0.val[k + 2N, l + 2N]*Gll.val[kPi + N, k]) + Gl0.val[k + 3N, l + 2N]*(2*Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
  Gl0.val[k, l + 2N]*Gll.val[kPi + N, k + 2N] + Gl0.val[k + N, l + 2N]*(2*Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) - 
  2*Gl0.val[kPi, l + 2N]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + Gl0.val[kPi + N, l + 2N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - 
    Gll.val[k + 3N, k + N]))*((G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l, kPi + 3N] - 
  (G00.val[l, lPj + 2N] + 2*G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi + 3N] + 2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*G0l.val[l + N, kPi + 3N] + 
  G00.val[l, lPj + 3N]*G0l.val[lPj + N, kPi + 3N] - (G00.val[l, lPj] + 2*G00.val[l + N, lPj + N])*G0l.val[lPj + 2N, kPi + 3N] + 
  G00.val[l, lPj + N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 
(-(Gl0.val[k + 2N, l + 3N]*Gll.val[kPi + N, k]) + Gl0.val[k + 3N, l + 3N]*(2*Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - Gl0.val[k, l + 3N]*Gll.val[kPi + N, k + 2N] + 
  Gl0.val[k + N, l + 3N]*(2*Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) - 2*Gl0.val[kPi, l + 3N]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + 
  Gl0.val[kPi + N, l + 3N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
 (2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*G0l.val[l, kPi + 3N] + G00.val[l + N, lPj + 2N]*G0l.val[lPj, kPi + 3N] + 
  (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*G0l.val[l + N, kPi + 3N] - 
  (2*G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi + 3N] + G00.val[l + N, lPj]*G0l.val[lPj + 2N, kPi + 3N] - 
  (2*G00.val[l, lPj] + G00.val[l + N, lPj + N])*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 
(2*Gl0.val[kPi + 3N, l + N]*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N]) + Gl0.val[k + 3N, l + N]*Gll.val[kPi + 2N, k + N] + Gl0.val[k + N, l + N]*Gll.val[kPi + 2N, k + 3N] + 
  Gl0.val[kPi + 2N, l + N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]) - 
  Gl0.val[k + 2N, l + N]*(Gll.val[kPi + 2N, k] + 2*Gll.val[kPi + 3N, k + N]) - Gl0.val[k, l + N]*(Gll.val[kPi + 2N, k + 2N] + 2*Gll.val[kPi + 3N, k + 3N]))*
 ((2*G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj + N, kPi] - 2*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj])*G0l.val[l + 2N, kPi] - 
  G00.val[l + 3N, lPj]*G0l.val[lPj + 2N, kPi] + (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l + 3N, kPi] + 
  (2*G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi] + G00.val[l + 3N, lPj + 2N]*(-G0l.val[lPj, kPi] + id*I[lPj, kPi])) - 
(-2*Gl0.val[kPi + 2N, l]*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) + Gl0.val[kPi + 3N, l]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]) - 
  Gl0.val[k + 2N, l]*Gll.val[kPi + 3N, k] + Gl0.val[k + 3N, l]*(2*Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - Gl0.val[k, l]*Gll.val[kPi + 3N, k + 2N] + 
  Gl0.val[k + N, l]*(2*Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*((G00.val[l + 2N, lPj + 2N] + 2*G00.val[l + 3N, lPj + 3N])*G0l.val[lPj, kPi + N] + 
  (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*G0l.val[l + 2N, kPi + N] + 
  (G00.val[l + 2N, lPj] + 2*G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi + N] - 2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*G0l.val[l + 3N, kPi + N] - 
  G00.val[l + 2N, lPj + N]*G0l.val[lPj + 3N, kPi + N] + G00.val[l + 2N, lPj + 3N]*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi]));
    
    return tcoscos1 + tcoscos2 +4tcoscos3 
end

@inline Base.@propagate_inbounds function full_nem_kernel(
    mc, model::TwoBandModel, klkPlP::NTuple{4}, packed_greens::_GM4{<: Matrix}, 
    flv, ::Union{Discrete_MBF1_X_symm, Discrete_MBF2_symm})
  k, l, kPi, lPj = klkPlP   #k, l, k+i, l+j
  G00, G0l, Gl0, Gll = packed_greens
  N = length(lattice(model))
  id = I[G0l.k, G0l.l] 

  tcoscos1=real(2*conj(Gll.val[kPi, k + N])*Gll.val[k, kPi + N] + 4*(conj(Gll.val[kPi + N, kPi]) + Gll.val[kPi, kPi + N])*(conj(Gll.val[k, k + N]) + Gll.val[k + N, k]) + 
  2*conj(Gll.val[kPi + N, k])*Gll.val[k + N, kPi] + 6*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N])*I[k, kPi] - 
  2*(Gll.val[kPi, k]*(Gll.val[k + N, kPi + N] + 2*real(Gll.val[k + N, kPi + N])) + Gll.val[k, kPi]*(Gll.val[kPi + N, k + N] + 2*real(Gll.val[kPi + N, k + N]))))*
 real(2*conj(G00.val[lPj, l + N])*G00.val[l, lPj + N] + 4*(conj(G00.val[lPj + N, lPj]) + G00.val[lPj, lPj + N])*(conj(G00.val[l, l + N]) + G00.val[l + N, l]) + 
  2*conj(G00.val[lPj + N, l])*G00.val[l + N, lPj] + 6*(G00.val[lPj, l] + G00.val[lPj + N, l + N])*I[l, lPj] - 
  2*(G00.val[lPj, l]*(G00.val[l + N, lPj + N] + 2*real(G00.val[l + N, lPj + N])) + G00.val[l, lPj]*(G00.val[lPj + N, l + N] + 2*real(G00.val[lPj + N, l + N]))));

  tcoscos2=(conj(Gl0.val[k, lPj + N])*conj(Gl0.val[kPi, l + N]) - conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi, lPj + N]))*
  (G0l.val[l, kPi + N]*G0l.val[lPj, k + N] - G0l.val[l, k + N]*G0l.val[lPj, kPi + N]) + 
 (conj(Gl0.val[k + N, lPj])*conj(Gl0.val[kPi + N, l]) - conj(Gl0.val[k + N, l])*conj(Gl0.val[kPi + N, lPj]))*
  (G0l.val[l + N, kPi]*G0l.val[lPj + N, k] - G0l.val[l + N, k]*G0l.val[lPj + N, kPi]) + 4*conj(G0l.val[l, k + N])*conj(Gl0.val[kPi + N, lPj])*G0l.val[lPj + N, kPi]*Gl0.val[k, l + N] + 
 4*(2*conj(G0l.val[lPj, kPi + N])*G0l.val[l, k + N] + conj(G0l.val[l, kPi + N])*G0l.val[lPj, k + N])*(2*conj(Gl0.val[kPi, lPj + N])*Gl0.val[k, l + N] + 
   conj(Gl0.val[k, lPj + N])*Gl0.val[kPi, l + N]) + (conj(G0l.val[l, kPi + N])*conj(G0l.val[lPj, k + N]) - conj(G0l.val[l, k + N])*conj(G0l.val[lPj, kPi + N]))*
  (Gl0.val[k, lPj + N]*Gl0.val[kPi, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi, lPj + N]) + 4*conj(G0l.val[lPj + N, kPi])*conj(Gl0.val[k, l + N])*G0l.val[l, k + N]*Gl0.val[kPi + N, lPj] + 
 4*(conj(G0l.val[lPj + N, k])*G0l.val[l + N, kPi] + 2*conj(G0l.val[l + N, k])*G0l.val[lPj + N, kPi])*(conj(Gl0.val[kPi + N, l])*Gl0.val[k + N, lPj] + 
   2*conj(Gl0.val[k + N, l])*Gl0.val[kPi + N, lPj]) + (conj(G0l.val[l + N, kPi])*conj(G0l.val[lPj + N, k]) - conj(G0l.val[l + N, k])*conj(G0l.val[lPj + N, kPi]))*
  (Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l] - Gl0.val[k + N, l]*Gl0.val[kPi + N, lPj]) - 
 4*G0l.val[lPj + N, kPi]*(conj(Gl0.val[kPi + N, lPj])*Gl0.val[k + N, l + N] + 2*conj(Gl0.val[k + N, lPj])*Gl0.val[kPi + N, l + N])*(conj(G0l.val[l, k]) - id*I[l, k]) - 
 4*conj(G0l.val[lPj + N, k])*(conj(Gl0.val[kPi + N, l + N])*Gl0.val[k + N, lPj] + 2*conj(Gl0.val[k + N, l + N])*Gl0.val[kPi + N, lPj])*(G0l.val[l, kPi] - id*I[l, kPi]) - 
 2*(Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + N])*(conj(G0l.val[lPj + N, kPi])*(conj(G0l.val[l, k]) - id*I[l, k]) + 
   conj(G0l.val[lPj + N, k])*(-conj(G0l.val[l, kPi]) + id*I[l, kPi])) - 
 2*(conj(Gl0.val[k + N, l + N])*conj(Gl0.val[kPi + N, lPj]) - conj(Gl0.val[k + N, lPj])*conj(Gl0.val[kPi + N, l + N]))*
  (G0l.val[lPj + N, kPi]*(G0l.val[l, k] - id*I[l, k]) + G0l.val[lPj + N, k]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
 4*conj(Gl0.val[kPi + N, l + N])*Gl0.val[k, lPj]*(G0l.val[l, kPi] - id*I[l, kPi])*(conj(G0l.val[lPj + N, k + N]) - id*I[lPj, k]) - 
 4*conj(G0l.val[l, kPi + N])*(2*conj(Gl0.val[kPi, lPj])*Gl0.val[k, l + N] + conj(Gl0.val[k, lPj])*Gl0.val[kPi, l + N])*(G0l.val[lPj + N, k + N] - id*I[lPj, k]) + 
 4*conj(Gl0.val[k, lPj])*Gl0.val[kPi + N, l + N]*(conj(G0l.val[l, kPi]) - id*I[l, kPi])*(G0l.val[lPj + N, k + N] - id*I[lPj, k]) - 
 4*conj(Gl0.val[k, l])*Gl0.val[kPi + N, lPj]*(conj(G0l.val[lPj + N, kPi])*(G0l.val[l + N, k + N] - id*I[l, k]) + 
   2*conj(G0l.val[l + N, kPi])*(G0l.val[lPj + N, k + N] - id*I[lPj, k])) - 2*(-(Gl0.val[k, lPj]*Gl0.val[kPi + N, l]) + Gl0.val[k, l]*Gl0.val[kPi + N, lPj])*
  (conj(G0l.val[lPj + N, kPi])*(conj(G0l.val[l + N, k + N]) - id*I[l, k]) + conj(G0l.val[l + N, kPi])*(-conj(G0l.val[lPj + N, k + N]) + id*I[lPj, k])) + 
 4*conj(Gl0.val[kPi + N, l])*Gl0.val[k, lPj]*(-2*G0l.val[lPj + N, kPi]*(conj(G0l.val[l + N, k + N]) - id*I[l, k]) + 
   G0l.val[l + N, kPi]*(-conj(G0l.val[lPj + N, k + N]) + id*I[lPj, k])) - 
 2*(-(conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi + N, l])) + conj(Gl0.val[k, l])*conj(Gl0.val[kPi + N, lPj]))*
  (G0l.val[lPj + N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[l + N, kPi]*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) + 
 4*(conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi + N, lPj]) - conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi + N, l + N]))*
  (4*conj(G0l.val[l, k + N])*conj(G0l.val[lPj + N, kPi]) + G0l.val[l, k + N]*G0l.val[lPj + N, kPi] + 4*(conj(G0l.val[l, kPi]) - id*I[l, kPi])*
    (-conj(G0l.val[lPj + N, k + N]) + id*I[lPj, k]) + (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) + 
 4*(Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + N, l + N])*(conj(G0l.val[l, k + N])*conj(G0l.val[lPj + N, kPi]) + 4*G0l.val[l, k + N]*G0l.val[lPj + N, kPi] + 
   (conj(G0l.val[l, kPi]) - id*I[l, kPi])*(-conj(G0l.val[lPj + N, k + N]) + id*I[lPj, k]) + 
   4*(G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) - 
 4*G0l.val[l, k + N]*(2*conj(Gl0.val[kPi, l + N])*Gl0.val[k, lPj] + conj(Gl0.val[k, l + N])*Gl0.val[kPi, lPj])*(conj(G0l.val[lPj + N, kPi + N]) - id*I[lPj, kPi]) + 
 2*(Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj + N])*(conj(G0l.val[lPj, k + N])*(-conj(G0l.val[l, kPi]) + id*I[l, kPi]) + 
   conj(G0l.val[l, k + N])*(conj(G0l.val[lPj, kPi]) - id*I[lPj, kPi])) + 4*conj(Gl0.val[k, lPj + N])*Gl0.val[kPi + N, l + N]*
  (G0l.val[lPj, k + N]*(-conj(G0l.val[l, kPi]) + id*I[l, kPi]) - 2*G0l.val[l, k + N]*(conj(G0l.val[lPj, kPi]) - id*I[lPj, kPi])) + 
 (-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj + N])*
  ((conj(G0l.val[l, kPi]) - id*I[l, kPi])*(-conj(G0l.val[lPj, k]) + id*I[lPj, k]) + (conj(G0l.val[l, k]) - id*I[l, k])*
    (conj(G0l.val[lPj, kPi]) - id*I[lPj, kPi])) + (-(Gl0.val[k, lPj]*Gl0.val[kPi, l]) + Gl0.val[k, l]*Gl0.val[kPi, lPj])*
  ((conj(G0l.val[l + N, kPi + N]) - id*I[l, kPi])*(-conj(G0l.val[lPj + N, k + N]) + id*I[lPj, k]) + 
   (conj(G0l.val[l + N, k + N]) - id*I[l, k])*(conj(G0l.val[lPj + N, kPi + N]) - id*I[lPj, kPi])) + 
 4*(2*conj(Gl0.val[kPi, l])*Gl0.val[k, lPj] + conj(Gl0.val[k, l])*Gl0.val[kPi, lPj])*
  (2*(conj(G0l.val[l + N, kPi + N]) - id*I[l, kPi])*(G0l.val[lPj + N, k + N] - id*I[lPj, k]) + 
   (G0l.val[l + N, k + N] - id*I[l, k])*(conj(G0l.val[lPj + N, kPi + N]) - id*I[lPj, kPi])) - 4*conj(Gl0.val[kPi + N, lPj + N])*Gl0.val[k, l + N]*
  (2*conj(G0l.val[lPj, k + N])*(G0l.val[l, kPi] - id*I[l, kPi]) + conj(G0l.val[l, k + N])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
 2*(conj(Gl0.val[k, lPj + N])*conj(Gl0.val[kPi + N, l + N]) - conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi + N, lPj + N]))*
  (G0l.val[lPj, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + G0l.val[l, k + N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
 4*(conj(Gl0.val[kPi + N, lPj + N])*Gl0.val[k + N, l + N] + 2*conj(Gl0.val[k + N, lPj + N])*Gl0.val[kPi + N, l + N])*
  (2*(G0l.val[l, kPi] - id*I[l, kPi])*(conj(G0l.val[lPj, k]) - id*I[lPj, k]) + (conj(G0l.val[l, k]) - id*I[l, k])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
 (-(conj(Gl0.val[k + N, lPj + N])*conj(Gl0.val[kPi + N, l + N])) + conj(Gl0.val[k + N, l + N])*conj(Gl0.val[kPi + N, lPj + N]))*
  ((G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj, k] + id*I[lPj, k]) + (G0l.val[l, k] - id*I[l, k])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
 (-(conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi, l])) + conj(Gl0.val[k, l])*conj(Gl0.val[kPi, lPj]))*
  ((G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k]) + (G0l.val[l + N, k + N] - id*I[l, k])*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
 2*(Gl0.val[k, l + N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + N])*(conj(G0l.val[l, kPi + N])*(conj(G0l.val[lPj + N, k + N]) - id*I[lPj, k]) + 
   conj(G0l.val[l, k + N])*(-conj(G0l.val[lPj + N, kPi + N]) + id*I[lPj, kPi])) + 
 2*(conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi, lPj]) - conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi, l + N]))*
  (G0l.val[l, kPi + N]*(G0l.val[lPj + N, k + N] - id*I[lPj, k]) + G0l.val[l, k + N]*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi]));

  tcoscos3=(-4*conj(Gl0.val[kPi, l + N])*(conj(Gll.val[k + N, k]) + Gll.val[k, k + N]) - 2*conj(Gl0.val[k, l + N])*Gll.val[kPi, k + N] + 
  conj(Gl0.val[k + N, l + N])*(4*conj(Gll.val[kPi, k]) + 2*Gll.val[kPi, k]))*(2*conj(G0l.val[lPj, kPi + N])*G00.val[l, lPj + N] + 
  4*conj(G0l.val[l, kPi + N])*(conj(G00.val[lPj + N, lPj]) + G00.val[lPj, lPj + N])  + 
  2*(conj(G0l.val[lPj + N, kPi + N]) - id*I[lPj, kPi])*(G00.val[l, lPj] - 4*real(G00.val[l, lPj]))) + 
(-4*conj(Gl0.val[kPi, l])*(conj(Gll.val[k + N, k]) + Gll.val[k, k + N]) - 2*conj(Gl0.val[k, l])*Gll.val[kPi, k + N] + 
  conj(Gl0.val[k + N, l])*(4*conj(Gll.val[kPi, k]) + 2*Gll.val[kPi, k]))*(2*conj(G0l.val[lPj + N, kPi + N])*G00.val[l + N, lPj] + 
  4*(conj(G00.val[lPj, lPj + N]) + G00.val[lPj + N, lPj])*(conj(G0l.val[l + N, kPi + N]) - id*I[l, kPi])  - 
  2*id*G00.val[l + N, lPj]*I[lPj, kPi] + 2*conj(G0l.val[lPj, kPi + N])*(G00.val[l + N, lPj + N] - 4*real(G00.val[l + N, lPj + N]))) + 
(-2*conj(G00.val[l + N, lPj])*G0l.val[lPj + N, kPi + N] - 4*(conj(G00.val[lPj + N, lPj]) + G00.val[lPj, lPj + N])*(G0l.val[l + N, kPi + N] - id*I[l, kPi])  
   + 2*id*conj(G00.val[l + N, lPj])*I[lPj, kPi] + 2*G0l.val[lPj, kPi + N]*(G00.val[l + N, lPj + N] + 2*real(G00.val[l + N, lPj + N])))*
 (2*conj(Gll.val[kPi, k + N])*Gl0.val[k, l] + 4*Gl0.val[kPi, l]*(conj(Gll.val[k, k + N]) + Gll.val[k + N, k]) + Gl0.val[k + N, l]*(-2*Gll.val[kPi, k]  - 4*real(Gll.val[kPi, k]))) - 
(4*(conj(G00.val[lPj, lPj + N]) + G00.val[lPj + N, lPj])*G0l.val[l, kPi + N] + 2*conj(G00.val[l, lPj + N])*G0l.val[lPj, kPi + N] 
- 2*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])*(G00.val[l, lPj] + 2*real(G00.val[l, lPj])))*
 (2*conj(Gll.val[kPi, k + N])*Gl0.val[k, l + N] + 4*Gl0.val[kPi, l + N]*(conj(Gll.val[k, k + N]) + Gll.val[k + N, k]) + 
  Gl0.val[k + N, l + N]*(-2*Gll.val[kPi, k]  - 4*real(Gll.val[kPi, k]))) - 
(2*conj(G0l.val[lPj + N, kPi])*G00.val[l + N, lPj] + 4*conj(G0l.val[l + N, kPi])*(conj(G00.val[lPj, lPj + N]) + G00.val[lPj + N, lPj]) 
+ 2*(conj(G0l.val[lPj, kPi]) - id*I[lPj, kPi])*(G00.val[l + N, lPj + N] - 4*real(G00.val[l + N, lPj + N])))*
 (4*conj(Gl0.val[kPi + N, l])*(conj(Gll.val[k, k + N]) + Gll.val[k + N, k]) + 2*conj(Gl0.val[k + N, l])*Gll.val[kPi + N, k] + 
  conj(Gl0.val[k, l])*(2*Gll.val[kPi + N, k + N]  - 8*real(Gll.val[kPi + N, k + N]))) - 
(2*conj(G0l.val[lPj, kPi])*G00.val[l, lPj + N] + 4*(conj(G00.val[lPj + N, lPj]) + G00.val[lPj, lPj + N])*(conj(G0l.val[l, kPi]) - id*I[l, kPi]) 
- 2*id*G00.val[l, lPj + N]*I[lPj, kPi] + 2*conj(G0l.val[lPj + N, kPi])*(G00.val[l, lPj] - 4*real(G00.val[l, lPj])))*
 (4*conj(Gl0.val[kPi + N, l + N])*(conj(Gll.val[k, k + N]) + Gll.val[k + N, k]) + 2*conj(Gl0.val[k + N, l + N])*Gll.val[kPi + N, k] + 
  conj(Gl0.val[k, l + N])*(2*Gll.val[kPi + N, k + N]  - 8*real(Gll.val[kPi + N, k + N]))) + 
(4*(conj(G00.val[lPj + N, lPj]) + G00.val[lPj, lPj + N])*G0l.val[l + N, kPi] + 2*conj(G00.val[l + N, lPj])*G0l.val[lPj + N, kPi]  - 
  2*(G0l.val[lPj, kPi] - id*I[lPj, kPi])*(G00.val[l + N, lPj + N] + 2*real(G00.val[l + N, lPj + N])))*(-2*conj(Gll.val[kPi + N, k])*Gl0.val[k + N, l] - 
  4*Gl0.val[kPi + N, l]*(conj(Gll.val[k + N, k]) + Gll.val[k, k + N]) + Gl0.val[k, l]*(2*Gll.val[kPi + N, k + N]  + 4*real(Gll.val[kPi + N, k + N]))) + 
(2*conj(G00.val[l, lPj + N])*G0l.val[lPj, kPi] + 4*(conj(G00.val[lPj, lPj + N]) + G00.val[lPj + N, lPj])*(G0l.val[l, kPi] - id*I[l, kPi])  - 
  2*id*conj(G00.val[l, lPj + N])*I[lPj, kPi] - 2*G0l.val[lPj + N, kPi]*(G00.val[l, lPj] + 2*real(G00.val[l, lPj])))*(-2*conj(Gll.val[kPi + N, k])*Gl0.val[k + N, l + N] - 
  4*Gl0.val[kPi + N, l + N]*(conj(Gll.val[k + N, k]) + Gll.val[k, k + N]) + Gl0.val[k, l + N]*(2*Gll.val[kPi + N, k + N]  + 4*real(Gll.val[kPi + N, k + N])));

  #Note that tcoscos3 was multiplied with 4 in kate
  return tcoscos1 + tcoscos2 +tcoscos3 
end




@inline Base.@propagate_inbounds function full_nem_kernel_old(mc::DQMC, model::TwoBandModel, klkPlP::NTuple{4}, 
  G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
  return full_nem_kernel_old(mc, model, klkPlP, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_nem_kernel_old(mc::DQMC, model::TwoBandModel, klkPlP::NTuple{4}, 
  G::GreensMatrix, flv, field::AbstractMagnBosonField)
  return full_nem_kernel_old(mc, model, klkPlP, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function full_nem_kernel_old(
  mc, model::TwoBandModel, klkPlP::NTuple{4}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
  )
  k, l, kPi, lPj = klkPlP   #k, l, k+i, l+j
  G00, G0l, Gl0, Gll = packed_greens
  N = length(lattice(model))
  id = I[G0l.k, G0l.l] 

  tcoscos1=(-(Gll.val[kPi, k + 3N]*Gll.val[k + N, kPi + 2N]) - Gll.val[k, kPi + 3N]*Gll.val[kPi + N, k + 2N] + Gll.val[k + N, kPi + 3N]*(2*Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
  Gll.val[k, kPi + 2N]*(Gll.val[kPi, k + 2N] + 2*Gll.val[kPi + N, k + 3N]) - Gll.val[kPi + N, k]*Gll.val[k + 2N, kPi + 3N] + 
  (Gll.val[k, kPi] + 2*Gll.val[k + N, kPi + N])*Gll.val[kPi + 2N, k + 2N] - Gll.val[k + N, kPi]*Gll.val[kPi + 2N, k + 3N] - 
  2*(Gll.val[kPi, kPi + 3N] + Gll.val[kPi + 2N, kPi + N])*(Gll.val[k + N, k + 2N] + Gll.val[k + 3N, k]) - Gll.val[kPi + 2N, k + N]*Gll.val[k + 3N, kPi] + 
  Gll.val[kPi + 2N, k]*(Gll.val[k + 2N, kPi] + 2*Gll.val[k + 3N, kPi + N]) - Gll.val[kPi, k + N]*Gll.val[k + 3N, kPi + 2N] - Gll.val[k + 2N, kPi + N]*Gll.val[kPi + 3N, k] - 
  2*(Gll.val[k, k + 3N] + Gll.val[k + 2N, k + N])*(Gll.val[kPi + N, kPi + 2N] + Gll.val[kPi + 3N, kPi]) + (2*Gll.val[k + 2N, kPi] + Gll.val[k + 3N, kPi + N])*Gll.val[kPi + 3N, k + N] - 
  (Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N])*(Gll.val[kPi, kPi + 2N] - Gll.val[kPi + N, kPi + 3N] + Gll.val[kPi + 2N, kPi] - 
  Gll.val[kPi + 3N, kPi + N]) - Gll.val[k, kPi + N]*Gll.val[kPi + 3N, k + 2N] + (2*Gll.val[k, kPi] + Gll.val[k + N, kPi + N])*Gll.val[kPi + 3N, k + 3N] + 
  Gll.val[kPi + N, k + N]*(2*Gll.val[k + 2N, kPi + 2N] + Gll.val[k + 3N, kPi + 3N] - 3*I[k, kPi]) + 
  Gll.val[kPi, k]*(Gll.val[k + 2N, kPi + 2N] + 2*Gll.val[k + 3N, kPi + 3N] - 3*I[k, kPi]) - 3*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N])*I[k, kPi])*
  (-(G00.val[lPj, l + 3N]*G00.val[l + N, lPj + 2N]) - G00.val[l, lPj + 3N]*G00.val[lPj + N, l + 2N] + G00.val[l + N, lPj + 3N]*(2*G00.val[lPj, l + 2N] + G00.val[lPj + N, l + 3N]) + 
  G00.val[l, lPj + 2N]*(G00.val[lPj, l + 2N] + 2*G00.val[lPj + N, l + 3N]) - G00.val[lPj + N, l]*G00.val[l + 2N, lPj + 3N] + 
  (G00.val[l, lPj] + 2*G00.val[l + N, lPj + N])*G00.val[lPj + 2N, l + 2N] - G00.val[l + N, lPj]*G00.val[lPj + 2N, l + 3N] - 
  2*(G00.val[lPj, lPj + 3N] + G00.val[lPj + 2N, lPj + N])*(G00.val[l + N, l + 2N] + G00.val[l + 3N, l]) - G00.val[lPj + 2N, l + N]*G00.val[l + 3N, lPj] + 
  G00.val[lPj + 2N, l]*(G00.val[l + 2N, lPj] + 2*G00.val[l + 3N, lPj + N]) - G00.val[lPj, l + N]*G00.val[l + 3N, lPj + 2N] - G00.val[l + 2N, lPj + N]*G00.val[lPj + 3N, l] - 
  2*(G00.val[l, l + 3N] + G00.val[l + 2N, l + N])*(G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 3N, lPj]) + (2*G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G00.val[lPj + 3N, l + N] - 
  (G00.val[l, l + 2N] - G00.val[l + N, l + 3N] + G00.val[l + 2N, l] - G00.val[l + 3N, l + N])*(G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - 
  G00.val[lPj + 3N, lPj + N]) - G00.val[l, lPj + N]*G00.val[lPj + 3N, l + 2N] + (2*G00.val[l, lPj] + G00.val[l + N, lPj + N])*G00.val[lPj + 3N, l + 3N] + 
  G00.val[lPj + N, l + N]*(2*G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N] - 3*I[l, lPj]) + 
  G00.val[lPj, l]*(G00.val[l + 2N, lPj + 2N] + 2*G00.val[l + 3N, lPj + 3N] - 3*I[l, lPj]) - 3*(G00.val[lPj + 2N, l + 2N] + G00.val[lPj + 3N, l + 3N])*I[l, lPj]);


  tcoscos2=(G0l.val[l + 3N, kPi + 2N]*G0l.val[lPj + 3N, k + 2N] - G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + 2N])*
  (Gl0.val[k, lPj + N]*Gl0.val[kPi, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi, lPj + N]) + 2*(G0l.val[l, kPi + 2N]*G0l.val[lPj + 3N, k + 2N] - G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + 2N])*
  (Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 2N]) + (G0l.val[l, kPi + 2N]*G0l.val[lPj, k + 2N] - G0l.val[l, k + 2N]*G0l.val[lPj, kPi + 2N])*
  (Gl0.val[k, lPj + 2N]*Gl0.val[kPi, l + 2N] - Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + 2N]) + 
  4*(-(G0l.val[l, kPi + 2N]*G0l.val[lPj + 3N, k + 2N]) + G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + 2N])*(Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + 3N]) + 
  2*(-(G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 3N, k + 2N]) + G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + 2N])*
  (Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 3N]) - 2*(G0l.val[lPj, kPi + 2N]*G0l.val[l + N, k + 2N] - G0l.val[lPj, k + 2N]*G0l.val[l + N, kPi + 2N])*
  (Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi, l + 3N]) + 
  4*(-(G0l.val[l, kPi + 2N]*G0l.val[lPj + N, k + 2N]) + G0l.val[l, k + 2N]*G0l.val[lPj + N, kPi + 2N])*(Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + 2N] - 
  Gl0.val[k, lPj + 2N]*Gl0.val[kPi, l + 3N]) + (G0l.val[l + N, kPi + 2N]*G0l.val[lPj + N, k + 2N] - G0l.val[l + N, k + 2N]*G0l.val[lPj + N, kPi + 2N])*
  (Gl0.val[k, lPj + 3N]*Gl0.val[kPi, l + 3N] - Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + 3N]) - 2*(G0l.val[l, k + 3N]*G0l.val[lPj, kPi + 2N] - G0l.val[l, kPi + 2N]*G0l.val[lPj, k + 3N])*
  (Gl0.val[kPi, lPj + 2N]*Gl0.val[k + N, l + 2N] - Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + 2N]) + 
  4*(G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi + 2N] - G0l.val[lPj, kPi + 2N]*G0l.val[l + N, k + 3N])*(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + 2N] - 
  Gl0.val[kPi, lPj + 2N]*Gl0.val[k + N, l + 3N]) + 8*(G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi + 2N] - G0l.val[l, kPi + 2N]*G0l.val[lPj + N, k + 3N])*
  (Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + 2N] - Gl0.val[kPi, lPj + 2N]*Gl0.val[k + N, l + 3N]) - 
  2*(G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi + 2N] - G0l.val[l + N, kPi + 2N]*G0l.val[lPj + N, k + 3N])*(Gl0.val[kPi, lPj + 3N]*Gl0.val[k + N, l + 3N] - 
  Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + 3N]) + (G0l.val[l + 2N, kPi + 3N]*G0l.val[lPj + 2N, k + 3N] - G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi + 3N])*
  (Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l] - Gl0.val[k + N, l]*Gl0.val[kPi + N, lPj]) + 2*(-(G0l.val[l, kPi + 3N]*G0l.val[lPj + 2N, k + 3N]) + G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi + 3N])*
  (Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 2N]) + 
  4*(-(G0l.val[l + N, kPi + 3N]*G0l.val[lPj + 2N, k + 3N]) + G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + 3N])*
  (Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 2N]) + 
  4*(G0l.val[l, k + 3N]*G0l.val[lPj, kPi + 2N] - G0l.val[l, kPi + 2N]*G0l.val[lPj, k + 3N])*(-(Gl0.val[k, lPj + 2N]*Gl0.val[kPi + N, l + 2N]) + 
  Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj + 2N]) + (G0l.val[l, kPi + 3N]*G0l.val[lPj, k + 3N] - G0l.val[l, k + 3N]*G0l.val[lPj, kPi + 3N])*
  (Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + N, l + 2N] - Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + 2N]) + 
  8*(G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi + 2N] - G0l.val[lPj, kPi + 2N]*G0l.val[l + N, k + 3N])*(Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj + 2N] - 
  Gl0.val[k, lPj + 2N]*Gl0.val[kPi + N, l + 3N]) + 16*(G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi + 2N] - G0l.val[l, kPi + 2N]*G0l.val[lPj + N, k + 3N])*
  (Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + N, l + 3N]) + 
  2*(G0l.val[l + N, kPi + 3N]*G0l.val[lPj + 2N, k + 3N] - G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + 3N])*
  (Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 3N]) - 
  2*(G0l.val[lPj, kPi + 3N]*G0l.val[l + N, k + 3N] - G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi + 3N])*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + 2N] - 
  Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + N, l + 3N]) + 4*(-(G0l.val[l, kPi + 3N]*G0l.val[lPj + N, k + 3N]) + G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi + 3N])*
  (Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + N, l + 3N]) + 
  4*(G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi + 2N] - G0l.val[l + N, kPi + 2N]*G0l.val[lPj + N, k + 3N])*(-(Gl0.val[k, lPj + 3N]*Gl0.val[kPi + N, l + 3N]) + 
  Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj + 3N]) + (G0l.val[l + N, kPi + 3N]*G0l.val[lPj + N, k + 3N] - G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi + 3N])*
  (Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + N, l + 3N] - Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + 3N]) + 
  2*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + 3N])*(Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l] - Gl0.val[k + N, l]*Gl0.val[kPi + 2N, lPj]) + 
  (G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k] - G0l.val[l + 2N, k]*G0l.val[lPj + 2N, kPi])*(Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l] - Gl0.val[k + 2N, l]*Gl0.val[kPi + 2N, lPj]) - 
  2*(G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k] - G0l.val[lPj + 2N, k]*G0l.val[l + 3N, kPi])*(Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + N]) + 
  4*(-(G0l.val[l + 2N, kPi]*G0l.val[lPj + 3N, k]) + G0l.val[l + 2N, k]*G0l.val[lPj + 3N, kPi])*(Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 2N, lPj] - 
  Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + N]) + 2*(G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k + 2N])*
  (-(Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + N]) + Gl0.val[k, l + N]*Gl0.val[kPi + 2N, lPj + N]) + 
  4*(G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + 2N])*(-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + N]) + 
  Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj + N]) + (G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k] - G0l.val[l + 3N, k]*G0l.val[lPj + 3N, kPi])*
  (Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + N] - Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 2N, lPj + N]) - 
  8*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + 3N])*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - 
  Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N]) - 8*(G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + 2N])*
  (Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N]) + 
  4*(-(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k]) + G0l.val[l + N, k]*G0l.val[lPj + 2N, kPi])*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - 
  Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + 2N]) + 4*(G0l.val[l, k + 2N]*G0l.val[lPj, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj, k + 2N])*
  (-(Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N]) + 
  4*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k + 2N])*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - 
  Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 3N]) + 4*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + 3N])*
  (Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N]) + 
  16*(G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + 2N])*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] - 
  Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N]) + 2*(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k] - G0l.val[l + N, k]*G0l.val[lPj + 2N, kPi])*
  (Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + 3N]) + 
  2*(-(G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k]) + G0l.val[l + N, k]*G0l.val[lPj + 3N, kPi])*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - 
  Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + 3N]) + 2*(G0l.val[l + N, k + 2N]*G0l.val[lPj + N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 2N])*
  (-(Gl0.val[k, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N]) + Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N]) + 
  2*(G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 3N])*(Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N] - 
  Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N]) + (G0l.val[l + N, kPi]*G0l.val[lPj + N, k] - G0l.val[l + N, k]*G0l.val[lPj + N, kPi])*
  (Gl0.val[k + 2N, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N]) - 
  2*(G0l.val[l + 2N, k + N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + N])*(Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l] - Gl0.val[kPi + 2N, l]*Gl0.val[k + 3N, lPj]) + 
  4*(G0l.val[lPj + 2N, k + N]*G0l.val[l + 3N, kPi] - G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k + N])*(Gl0.val[kPi + 2N, l + N]*Gl0.val[k + 3N, lPj] - 
  Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + N]) + 8*(G0l.val[l + 2N, k + N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 3N, k + N])*
  (Gl0.val[kPi + 2N, l + N]*Gl0.val[k + 3N, lPj] - Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + N]) - 
  2*(G0l.val[l + 3N, k + N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k + N])*(Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + N] - 
  Gl0.val[kPi + 2N, l + N]*Gl0.val[k + 3N, lPj + N]) + 4*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + 3N])*
  (-(Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l]) + Gl0.val[k, l]*Gl0.val[kPi + 3N, lPj]) + 
  2*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 2N, k + 3N])*
  (-(Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l]) + Gl0.val[k + N, l]*Gl0.val[kPi + 3N, lPj]) + 4*(G0l.val[l + 2N, k + N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + N])*
  (-(Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l]) + Gl0.val[k + 2N, l]*Gl0.val[kPi + 3N, lPj]) + 
  (G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 2N, k + N] - G0l.val[l + 2N, k + N]*G0l.val[lPj + 2N, kPi + N])*(Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l] - 
  Gl0.val[k + 3N, l]*Gl0.val[kPi + 3N, lPj]) + 8*(G0l.val[lPj + 2N, k + N]*G0l.val[l + 3N, kPi] - G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k + N])*
  (Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + N]) + 
  16*(G0l.val[l + 2N, k + N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 3N, k + N])*(Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 3N, lPj] - 
  Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + N]) - 2*(G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + N] - G0l.val[lPj + 2N, k + N]*G0l.val[l + 3N, kPi + N])*
  (Gl0.val[k + 3N, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + N]) + 
  4*(-(G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 3N, k + N]) + G0l.val[l + 2N, k + N]*G0l.val[lPj + 3N, kPi + N])*
  (Gl0.val[k + 3N, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + N]) + 
  2*(G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + 2N])*
  (Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj + N]) + 
  4*(G0l.val[l + 3N, k + N]*G0l.val[lPj + 3N, kPi] - G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k + N])*(-(Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + N]) + 
  Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 3N, lPj + N]) + (G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + N] - G0l.val[l + 3N, k + N]*G0l.val[lPj + 3N, kPi + N])*
  (Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + N] - Gl0.val[k + 3N, l + N]*Gl0.val[kPi + 3N, lPj + N]) + 
  16*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + 3N])*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - 
  Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N]) + 4*(G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + 2N])*
  (Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N]) + 
  4*(G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj + 2N, k + 3N])*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj] - 
  Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + 2N]) + 2*(-(G0l.val[l, kPi + N]*G0l.val[lPj + 2N, k + N]) + G0l.val[l, k + N]*G0l.val[lPj + 2N, kPi + N])*
  (Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + 2N]) + 
  2*(G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + N] - G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi + N])*(Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - 
  Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + 2N]) + 2*(G0l.val[l, k + 2N]*G0l.val[lPj, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj, k + 2N])*
  (Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N] - Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N]) + 2*(G0l.val[l, k + 3N]*G0l.val[lPj, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj, k + 3N])*
  (-(Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N]) + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N]) + 
  (G0l.val[l, kPi + N]*G0l.val[lPj, k + N] - G0l.val[l, k + N]*G0l.val[lPj, kPi + N])*(Gl0.val[k + 3N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N] - 
  Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N]) - 8*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + 3N])*
  (Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N]) - 8*(G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + 2N])*
  (Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N]) + 4*(-(G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + N]) + G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi + N])*
  (Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + 3N]) + 
  4*(G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi] - G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 3N])*(-(Gl0.val[k, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + 
  Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N]) - 8*(-(Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + N]) + Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 2N])*
  (-(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + N]) + G0l.val[lPj + 2N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k])) + 
  4*(-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj]) + Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 3N])*(-(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + N]) + 
  G0l.val[lPj + 2N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k])) + 16*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*
  (-(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + N]) + G0l.val[lPj + 2N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k])) - 
  4*(-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + N]) + Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 3N])*
  (-(G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k + N]) + G0l.val[lPj + 3N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k])) + 
  8*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*
  (-(G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k + N]) + G0l.val[lPj + 3N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k])) + 
  8*(Gl0.val[k, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + N])*(-(G0l.val[l + 2N, kPi]*G0l.val[lPj + 3N, k + 2N]) + 
  G0l.val[lPj + 3N, kPi]*(G0l.val[l + 2N, k + 2N] - id*I[l, k])) + 16*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + N])*
  (-(G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 3N, k + 2N]) + G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k])) - 
  8*(Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + N])*(-(G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 3N, k + 2N]) + 
  G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k])) + 4*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + N])*
  (-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi]) + G0l.val[lPj + 2N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k])) - 
  4*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + N])*(-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi + N]) + 
  G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k])) + 8*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + N] + G0l.val[lPj + 2N, kPi]*(-G0l.val[l + N, k + N] + id*I[l, k])) + 
  8*(Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + N])*(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi] + 
  G0l.val[lPj + 2N, kPi]*(-G0l.val[l + 3N, k + 3N] + id*I[l, k])) + 4*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
  (G0l.val[lPj + N, kPi]*(G0l.val[l, k] - id*I[l, k]) + G0l.val[lPj + N, k]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
  8*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*(G0l.val[l, k + 2N]*G0l.val[lPj + N, kPi] + 
  G0l.val[lPj + N, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 8*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi] + G0l.val[lPj + N, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
  16*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + N, kPi] + G0l.val[lPj + N, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
  2*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + 2N])*(G0l.val[lPj + 2N, kPi]*(G0l.val[l, k] - id*I[l, k]) + 
  G0l.val[lPj + 2N, k]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 4*(-(Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj]) + Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 2N])*
  (G0l.val[l, k + N]*G0l.val[lPj + 2N, kPi] + G0l.val[lPj + 2N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
  8*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + 2N])*(G0l.val[l, k + N]*G0l.val[lPj + 2N, kPi] + 
  G0l.val[lPj + 2N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 4*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 2N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi] + G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
  8*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 2N])*(G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi] + 
  G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 2*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + 2N])*
  (G0l.val[lPj + 3N, kPi]*(G0l.val[l, k] - id*I[l, k]) + G0l.val[lPj + 3N, k]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
  4*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + 3N])*(G0l.val[lPj + 3N, kPi]*(G0l.val[l, k] - id*I[l, k]) + 
  G0l.val[lPj + 3N, k]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 4*(-(Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + N]) + Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 2N])*
  (G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi] + G0l.val[lPj + 3N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 
  8*(-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj]) + Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 3N])*
  (G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi] + G0l.val[lPj + 3N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 
  8*(Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*
  (G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi] + G0l.val[lPj + 3N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 
  16*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi] + G0l.val[lPj + 3N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 
  4*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 2N])*(G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi] + 
  G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + 8*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 3N])*
  (G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi] + G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi])) - 
  2*(Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + 3N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[lPj, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[lPj, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) - 
  8*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
  (G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 2N] + G0l.val[lPj, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 
  4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 2N] + 
  G0l.val[lPj, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) - 4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 3N] + G0l.val[lPj, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 
  4*(Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*(G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + 
  G0l.val[lPj + 2N, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) - 2*(Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[lPj + 2N, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 
  8*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + N] + 
  G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) - 4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + N] + G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 
  2*(Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + 
  G0l.val[lPj + 3N, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 8*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 3N])*
  (G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) - 
  4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + N] + 
  G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + 4*(Gl0.val[k, l + N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + N])*
  (G0l.val[lPj + 3N, kPi + 2N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l + 2N, kPi + 2N] + id*I[l, kPi])) - 
  2*(Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + N])*(G0l.val[lPj + 2N, kPi + 3N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + 
  G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + 3N, kPi + 3N] + id*I[l, kPi])) - 
  2*(Gl0.val[kPi + 2N, lPj + 3N]*Gl0.val[k + 3N, l + 3N] - Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + 3N])*
  (G0l.val[lPj + N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[l + N, kPi]*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) + 
  4*(-(Gl0.val[k + 2N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*
  (G0l.val[lPj + N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[l + N, kPi]*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) - 
  8*(-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + 2N]) + Gl0.val[kPi + 2N, lPj + 2N]*Gl0.val[k + 3N, l + 3N])*
  (G0l.val[l, k + N]*G0l.val[lPj + N, kPi] + (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) + 
  16*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l, k + N]*G0l.val[lPj + N, kPi] + (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k])) + 
  8*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 2N])*(G0l.val[l, k + 2N]*G0l.val[lPj + 2N, kPi + N] + 
  G0l.val[l, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 4*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 2N])*
  (G0l.val[l, k + 2N]*G0l.val[lPj + 2N, kPi + N] + G0l.val[l, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 
  8*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 2N])*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi] + 
  G0l.val[l + N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 3N])*
  (G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi] + G0l.val[l + N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 
  2*(-(Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l]) + Gl0.val[k, l]*Gl0.val[kPi + 2N, lPj])*(G0l.val[lPj + 2N, kPi]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
  G0l.val[l + 2N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 4*(-(Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l]) + Gl0.val[k + N, l]*Gl0.val[kPi + 2N, lPj])*
  (G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + G0l.val[l + 2N, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 
  2*(-(Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l]) + Gl0.val[k, l]*Gl0.val[kPi + 3N, lPj])*(G0l.val[lPj + 2N, kPi + N]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
  G0l.val[l + 2N, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 4*(Gl0.val[k, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + N])*
  (G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k + 2N] + G0l.val[l + 3N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 
  8*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + N])*(G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + 2N] + 
  G0l.val[l + 3N, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 4*(Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + N])*
  (G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + 2N] + G0l.val[l + 3N, kPi + N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 
  4*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 2N])*(G0l.val[l, k + 2N]*G0l.val[lPj + 2N, kPi] + 
  (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 
  16*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N])*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi + N] + 
  (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 
  8*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N])*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi + N] + 
  (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 
  8*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi + N] + 
  (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) + 4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi + N] + (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 
  4*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*(G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + N] + 
  G0l.val[l, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 8*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[l, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  4*(-(Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + N]) + Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 2N])*(G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + 
  G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 8*(-(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 3N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  8*(Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + 2N])*(G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + 
  G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 16*(Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + N, l + 3N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi] + 
  G0l.val[l + N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 8*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi] + G0l.val[l + N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  4*(-(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + N]) + Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + 
  G0l.val[l + N, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 8*(Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + 3N])*
  (G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + G0l.val[l + N, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  8*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + N])*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi] + 
  G0l.val[l + 2N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 16*(Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + N])*
  (G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi] + G0l.val[l + 2N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  8*(Gl0.val[k + N, l + N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + N])*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi + N] + 
  G0l.val[l + 2N, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 2*(-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj + N])*
  (G0l.val[lPj + 3N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + G0l.val[l + 3N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  4*(-(Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + N]) + Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj + N])*(G0l.val[lPj + 3N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + 
  G0l.val[l + 3N, kPi]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 2*(-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + 3N, lPj + N])*
  (G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + G0l.val[l + 3N, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  2*(Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + N] - Gl0.val[kPi, l + N]*Gl0.val[k + N, lPj + N])*(G0l.val[lPj + 3N, kPi + 2N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + 
  G0l.val[l + 3N, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 4*(-(Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + N]) + Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj + N])*
  (G0l.val[lPj + 3N, kPi + 2N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + G0l.val[l + 3N, kPi + 2N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  4*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi] + (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  8*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N])*(G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi] + 
  (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 8*(Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N])*
  (G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi] + (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  16*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi] + 
  (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi + N] + 
  (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 8*(-(Gl0.val[kPi, l + N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + N])*
  (G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + (G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  16*(Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + N, l + N])*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + 
  (G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  2*(Gl0.val[kPi + 2N, lPj + 2N]*Gl0.val[k + 3N, l + 2N] - Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + 2N])*
  (G0l.val[lPj, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + G0l.val[l, k + N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
  4*(-(Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N]) + Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N])*
  (G0l.val[lPj, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + G0l.val[l, k + N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
  2*(-(Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N])*(G0l.val[lPj, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
  G0l.val[l, k + 2N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 2*(-(Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N])*
  (G0l.val[lPj, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + G0l.val[l, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
  4*(-(Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N]) + Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N])*(G0l.val[lPj, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
  G0l.val[l, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 2*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
  (G0l.val[l + N, kPi]*(-G0l.val[lPj, k] + id*I[lPj, k]) + G0l.val[l + N, k]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
  4*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*(-(G0l.val[lPj, k + 2N]*G0l.val[l + N, kPi]) + 
  G0l.val[l + N, k + 2N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
  (-(G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
  8*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*(-(G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi]) + 
  G0l.val[l + N, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + (-(Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N])*
  ((G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj, k] + id*I[lPj, k]) + (G0l.val[l, k] - id*I[l, k])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
  4*(-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + 2N]) + Gl0.val[kPi + 2N, lPj + 2N]*Gl0.val[k + 3N, l + 3N])*
  (-(G0l.val[lPj, k + N]*G0l.val[l + N, kPi]) + (G0l.val[l + N, k + N] - id*I[l, k])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
  8*(Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (-(G0l.val[lPj, k + N]*G0l.val[l + N, kPi]) + (G0l.val[l + N, k + N] - id*I[l, k])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
  4*(Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + 3N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l, kPi + N]*(-G0l.val[lPj + N, k + N] + id*I[lPj, k]) + G0l.val[l, k + N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
  16*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N])*
  (-(G0l.val[l, kPi + N]*G0l.val[lPj + N, k + 2N]) + G0l.val[l, k + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
  8*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (-(G0l.val[l, kPi + N]*G0l.val[lPj + N, k + 3N]) + G0l.val[l, k + 3N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
  4*(-(Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N]) + Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N])*(G0l.val[lPj + N, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + 
  G0l.val[l + N, k + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) - 2*(-(Gl0.val[k, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*
  (G0l.val[lPj + N, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + G0l.val[l + N, k + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
  2*(-(Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*(G0l.val[lPj + N, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + 
  G0l.val[l + N, k + 3N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
  (-(Gl0.val[k + 3N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N]) + Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*
  ((G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k]) + (G0l.val[l + N, k + N] - id*I[l, k])*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
  2*(Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + 2N])*(G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
  G0l.val[l, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 4*(-(Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 2N])*
  (-(G0l.val[l, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) + G0l.val[l, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
  8*(Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + N, l + 2N])*(-(G0l.val[l, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) + 
  G0l.val[l, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 4*(Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 2N])*
  (G0l.val[l + N, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l + N, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
  2*(Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + 3N])*(G0l.val[l + N, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
  G0l.val[l + N, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 8*(-(Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + N]) + Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 2N])*
  (-(G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
  4*(-(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 3N])*(-(G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) + 
  G0l.val[l + N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 16*(Gl0.val[k, l + 2N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + 2N])*
  (-(G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
  8*(Gl0.val[k, l + 3N]*Gl0.val[kPi + N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + N, l + 3N])*(-(G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) + 
  G0l.val[l + N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 2*(Gl0.val[kPi, lPj]*Gl0.val[k + N, l] - Gl0.val[kPi, l]*Gl0.val[k + N, lPj])*
  (G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + 2N, kPi + 2N] + id*I[l, kPi]) + G0l.val[l + 2N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
  4*(-(Gl0.val[k, lPj]*Gl0.val[kPi + N, l]) + Gl0.val[k, l]*Gl0.val[kPi + N, lPj])*(G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + 2N, kPi + 2N] + id*I[l, kPi]) + 
  G0l.val[l + 2N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 2*(Gl0.val[k, l + N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + N])*
  (G0l.val[l + 3N, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l + 3N, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
  (-(Gl0.val[k, lPj]*Gl0.val[kPi, l]) + Gl0.val[k, l]*Gl0.val[kPi, lPj])*((G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
  (G0l.val[l + 2N, k + 2N] - id*I[l, k])*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 4*(-(Gl0.val[kPi, l + N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + N])*
  (-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi + 2N]) + (G0l.val[l + 3N, k + 3N] - id*I[l, k])*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
  8*(Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + N, l + N])*(-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi + 2N]) + 
  (G0l.val[l + 3N, k + 3N] - id*I[l, k])*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
  2*(Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 2N])*(G0l.val[l, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + 
  G0l.val[l, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 4*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 3N])*
  (G0l.val[l, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + G0l.val[l, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 
  2*(Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 3N])*(G0l.val[l + N, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + 
  G0l.val[l + N, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 4*(Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + N])*
  (G0l.val[l + 2N, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + G0l.val[l + 2N, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 
  (-(Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj + N])*
  ((G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + (G0l.val[l + 3N, k + 3N] - id*I[l, k])*
  (G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 8*(Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*
  (G0l.val[l, kPi + N]*G0l.val[lPj + N, k + 2N] + G0l.val[l, k + 2N]*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi]));

  tcoscos3=(G00.val[l, l + 3N] - G00.val[l + N, l + 2N] + G00.val[l + 2N, l + N] - G00.val[l + 3N, l])*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N])*
  (G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi, lPj] - G0l.val[lPj + 2N, kPi + 2N]*Gl0.val[kPi, lPj + N] + G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi, lPj + 2N] - 
  G0l.val[lPj, kPi + 2N]*Gl0.val[kPi, lPj + 3N] - G0l.val[lPj + 3N, kPi + 3N]*Gl0.val[kPi + N, lPj] + G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi + N, lPj + N] - 
  G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi + N, lPj + 2N] + G0l.val[lPj, kPi + 3N]*Gl0.val[kPi + N, lPj + 3N] + G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 2N, lPj] - 
  G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 2N, lPj + N] + G0l.val[lPj + N, kPi]*Gl0.val[kPi + 2N, lPj + 2N] - G0l.val[lPj, kPi]*Gl0.val[kPi + 2N, lPj + 3N] - 
  G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 3N, lPj] + G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 3N, lPj + N] - G0l.val[lPj + N, kPi + N]*Gl0.val[kPi + 3N, lPj + 2N] + 
  G0l.val[lPj, kPi + N]*Gl0.val[kPi + 3N, lPj + 3N] + id*(Gl0.val[kPi, lPj + N] + Gl0.val[kPi + N, lPj] + Gl0.val[kPi + 2N, lPj + 3N] + Gl0.val[kPi + 3N, lPj + 2N])*I[lPj, kPi]) - 
  (-G00.val[l, l + 2N] + G00.val[l + N, l + 3N] - G00.val[l + 2N, l] + G00.val[l + 3N, l + N])*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k])*
  (G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi, lPj] - G0l.val[lPj + 3N, kPi + 3N]*Gl0.val[kPi, lPj + N] + G0l.val[lPj, kPi + 3N]*Gl0.val[kPi, lPj + 2N] - 
  G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi, lPj + 3N] - G0l.val[lPj + 2N, kPi + 2N]*Gl0.val[kPi + N, lPj] + G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi + N, lPj + N] - 
  G0l.val[lPj, kPi + 2N]*Gl0.val[kPi + N, lPj + 2N] + G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi + N, lPj + 3N] + G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 2N, lPj] - 
  G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 2N, lPj + N] + G0l.val[lPj, kPi + N]*Gl0.val[kPi + 2N, lPj + 2N] - G0l.val[lPj + N, kPi + N]*Gl0.val[kPi + 2N, lPj + 3N] - 
  G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 3N, lPj] + G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 3N, lPj + N] - G0l.val[lPj, kPi]*Gl0.val[kPi + 3N, lPj + 2N] + 
  G0l.val[lPj + N, kPi]*Gl0.val[kPi + 3N, lPj + 3N] + id*(Gl0.val[kPi, lPj + N] + Gl0.val[kPi + N, lPj] + Gl0.val[kPi + 2N, lPj + 3N] + Gl0.val[kPi + 3N, lPj + 2N])*I[lPj, kPi]) - 
  (G00.val[l, l + 3N] + G00.val[l + N, l + 2N] + G00.val[l + 2N, l + N] + G00.val[l + 3N, l])*(Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k])*
  (G0l.val[lPj + 3N, kPi + 3N]*Gl0.val[kPi, lPj] + G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi, lPj + N] + G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi, lPj + 2N] + 
  G0l.val[lPj, kPi + 3N]*Gl0.val[kPi, lPj + 3N] + G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi + N, lPj] + G0l.val[lPj + 2N, kPi + 2N]*Gl0.val[kPi + N, lPj + N] + 
  G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi + N, lPj + 2N] + G0l.val[lPj, kPi + 2N]*Gl0.val[kPi + N, lPj + 3N] + G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 2N, lPj] + 
  G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 2N, lPj + N] + G0l.val[lPj + N, kPi + N]*Gl0.val[kPi + 2N, lPj + 2N] + G0l.val[lPj, kPi + N]*Gl0.val[kPi + 2N, lPj + 3N] + 
  G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 3N, lPj] + G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 3N, lPj + N] + G0l.val[lPj + N, kPi]*Gl0.val[kPi + 3N, lPj + 2N] + 
  G0l.val[lPj, kPi]*Gl0.val[kPi + 3N, lPj + 3N] - id*Gl0.val[kPi, lPj]*I[lPj, kPi] - id*Gl0.val[kPi + N, lPj + N]*I[lPj, kPi] - id*Gl0.val[kPi + 2N, lPj + 2N]*I[lPj, kPi] - 
  id*Gl0.val[kPi + 3N, lPj + 3N]*I[lPj, kPi]) + (G00.val[l, l + 3N] + G00.val[l + N, l + 2N] + G00.val[l + 2N, l + N] + G00.val[l + 3N, l])*
  (Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k])*(G0l.val[lPj + 3N, kPi + 3N]*Gl0.val[kPi, lPj] + 
  G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi, lPj + N] + G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi, lPj + 2N] + G0l.val[lPj, kPi + 3N]*Gl0.val[kPi, lPj + 3N] - 
  G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi + N, lPj] - G0l.val[lPj + 2N, kPi + 2N]*Gl0.val[kPi + N, lPj + N] - G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi + N, lPj + 2N] - 
  G0l.val[lPj, kPi + 2N]*Gl0.val[kPi + N, lPj + 3N] + G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 2N, lPj] + G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 2N, lPj + N] + 
  G0l.val[lPj + N, kPi + N]*Gl0.val[kPi + 2N, lPj + 2N] + G0l.val[lPj, kPi + N]*Gl0.val[kPi + 2N, lPj + 3N] - G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 3N, lPj] - 
  G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 3N, lPj + N] - G0l.val[lPj + N, kPi]*Gl0.val[kPi + 3N, lPj + 2N] - G0l.val[lPj, kPi]*Gl0.val[kPi + 3N, lPj + 3N] - 
  id*Gl0.val[kPi, lPj]*I[lPj, kPi] + id*Gl0.val[kPi + N, lPj + N]*I[lPj, kPi] - id*Gl0.val[kPi + 2N, lPj + 2N]*I[lPj, kPi] + id*Gl0.val[kPi + 3N, lPj + 3N]*I[lPj, kPi]) + 
  (G00.val[l, l + 3N] - G00.val[l + N, l + 2N] + G00.val[l + 2N, l + N] - G00.val[l + 3N, l])*(-Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k])*
  (G0l.val[lPj + 3N, kPi + 3N]*Gl0.val[kPi, lPj] - G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi, lPj + N] + G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi, lPj + 2N] - 
  G0l.val[lPj, kPi + 3N]*Gl0.val[kPi, lPj + 3N] + G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi + N, lPj] - G0l.val[lPj + 2N, kPi + 2N]*Gl0.val[kPi + N, lPj + N] + 
  G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi + N, lPj + 2N] - G0l.val[lPj, kPi + 2N]*Gl0.val[kPi + N, lPj + 3N] + G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 2N, lPj] - 
  G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 2N, lPj + N] + G0l.val[lPj + N, kPi + N]*Gl0.val[kPi + 2N, lPj + 2N] - G0l.val[lPj, kPi + N]*Gl0.val[kPi + 2N, lPj + 3N] + 
  G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 3N, lPj] - G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 3N, lPj + N] + G0l.val[lPj + N, kPi]*Gl0.val[kPi + 3N, lPj + 2N] - 
  G0l.val[lPj, kPi]*Gl0.val[kPi + 3N, lPj + 3N] + id*(-Gl0.val[kPi, lPj] + Gl0.val[kPi + N, lPj + N] - Gl0.val[kPi + 2N, lPj + 2N] + Gl0.val[kPi + 3N, lPj + 3N])*I[lPj, kPi]) + 
  (-G00.val[l, l + 3N] - G00.val[l + N, l + 2N] - G00.val[l + 2N, l + N] - G00.val[l + 3N, l])*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N])*
  (-(G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi, lPj]) - G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi, lPj + 2N] - G0l.val[lPj, kPi + 2N]*Gl0.val[kPi, lPj + 3N] + 
  G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi + N, lPj + N] + G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi + N, lPj + 2N] + G0l.val[lPj, kPi + 3N]*Gl0.val[kPi + N, lPj + 3N] - 
  G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 2N, lPj] - G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 2N, lPj + N] - G0l.val[lPj + N, kPi]*Gl0.val[kPi + 2N, lPj + 2N] + 
  G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 3N, lPj] + G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 3N, lPj + N] + G0l.val[lPj, kPi + N]*Gl0.val[kPi + 3N, lPj + 3N] + 
  Gl0.val[kPi + 3N, lPj + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi]) + Gl0.val[kPi + N, lPj]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi]) + 
  Gl0.val[kPi + 2N, lPj + 3N]*(-G0l.val[lPj, kPi] + id*I[lPj, kPi]) + Gl0.val[kPi, lPj + N]*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi])) + 
  (-G00.val[l, l + 2N] + G00.val[l + N, l + 3N] - G00.val[l + 2N, l] + G00.val[l + 3N, l + N])*(-Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k])*
  (-(G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi, lPj]) - G0l.val[lPj, kPi + 3N]*Gl0.val[kPi, lPj + 2N] + G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi, lPj + 3N] + 
  G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi + N, lPj + N] - G0l.val[lPj, kPi + 2N]*Gl0.val[kPi + N, lPj + 2N] + G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi + N, lPj + 3N] - 
  G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 2N, lPj] + G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 2N, lPj + N] - G0l.val[lPj, kPi + N]*Gl0.val[kPi + 2N, lPj + 2N] - 
  G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 3N, lPj] + G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 3N, lPj + N] + G0l.val[lPj + N, kPi]*Gl0.val[kPi + 3N, lPj + 3N] + 
  Gl0.val[kPi + 2N, lPj + 3N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi]) + Gl0.val[kPi, lPj + N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi]) + 
  Gl0.val[kPi + 3N, lPj + 2N]*(-G0l.val[lPj, kPi] + id*I[lPj, kPi]) + Gl0.val[kPi + N, lPj]*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi])) - 
  (G00.val[l, l + 3N] - G00.val[l + N, l + 2N] + G00.val[l + 2N, l + N] - G00.val[l + 3N, l])*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k])*
  (G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi, lPj + N] - G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi, lPj + 2N] + G0l.val[lPj, kPi + 3N]*Gl0.val[kPi, lPj + 3N] + 
  G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi + N, lPj] + G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi + N, lPj + 2N] - G0l.val[lPj, kPi + 2N]*Gl0.val[kPi + N, lPj + 3N] - 
  G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 2N, lPj] + G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 2N, lPj + N] + G0l.val[lPj, kPi + N]*Gl0.val[kPi + 2N, lPj + 3N] + 
  G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 3N, lPj] - G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 3N, lPj + N] + G0l.val[lPj + N, kPi]*Gl0.val[kPi + 3N, lPj + 2N] + 
  Gl0.val[kPi + 3N, lPj + 3N]*(-G0l.val[lPj, kPi] + id*I[lPj, kPi]) + Gl0.val[kPi + 2N, lPj + 2N]*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi]) + 
  Gl0.val[kPi + N, lPj + N]*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi]) + Gl0.val[kPi, lPj]*(-G0l.val[lPj + 3N, kPi + 3N] + id*I[lPj, kPi])) + 
  (-G00.val[l, l + 2N] + G00.val[l + N, l + 3N] - G00.val[l + 2N, l] + G00.val[l + 3N, l + N])*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N])*
  (G0l.val[lPj + 3N, kPi + 2N]*Gl0.val[kPi, lPj + N] - G0l.val[lPj, kPi + 2N]*Gl0.val[kPi, lPj + 2N] + G0l.val[lPj + N, kPi + 2N]*Gl0.val[kPi, lPj + 3N] + 
  G0l.val[lPj + 2N, kPi + 3N]*Gl0.val[kPi + N, lPj] + G0l.val[lPj, kPi + 3N]*Gl0.val[kPi + N, lPj + 2N] - G0l.val[lPj + N, kPi + 3N]*Gl0.val[kPi + N, lPj + 3N] - 
  G0l.val[lPj + 2N, kPi]*Gl0.val[kPi + 2N, lPj] + G0l.val[lPj + 3N, kPi]*Gl0.val[kPi + 2N, lPj + N] + G0l.val[lPj + N, kPi]*Gl0.val[kPi + 2N, lPj + 3N] + 
  G0l.val[lPj + 2N, kPi + N]*Gl0.val[kPi + 3N, lPj] - G0l.val[lPj + 3N, kPi + N]*Gl0.val[kPi + 3N, lPj + N] + G0l.val[lPj, kPi + N]*Gl0.val[kPi + 3N, lPj + 2N] + 
  Gl0.val[kPi + 2N, lPj + 2N]*(-G0l.val[lPj, kPi] + id*I[lPj, kPi]) + Gl0.val[kPi + 3N, lPj + 3N]*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi]) + 
  Gl0.val[kPi, lPj]*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi]) + Gl0.val[kPi + N, lPj + N]*(-G0l.val[lPj + 3N, kPi + 3N] + id*I[lPj, kPi]));


  return tcoscos1 + tcoscos2 +4tcoscos3 
end


