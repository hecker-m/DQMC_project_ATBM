

###########################
### heat capacity h₂ measurement 
###########################
function heat_cap_h2(
    dqmc::DQMC, model::Model;
        greens_iterator = Greens(),
        lattice_iterator = EachDistancedBondPairSummed(dqmc), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = heat_cap_h2_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end






###########################
### heat capacity h₃ measurement 
###########################
function heat_cap_h3(
    dqmc::DQMC, model::Model;
        greens_iterator = Greens(),
        lattice_iterator = EachBondEachSiteSummed(dqmc), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = heat_cap_h3_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end





###########################
### heat capacity h₄ measurement 
###########################
function heat_cap_h4(
    dqmc::DQMC, model::Model;
        greens_iterator = Greens(),
        lattice_iterator = EachSiteTwiceSummed(), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = heat_cap_h4_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end












###########################
### heat capacity kernels for h₂, h₃, and h₄
###########################



###########################
### heat capacity kernel for h₂
###########################


"""
Calculates the heat capacity kernel for h₂
"""
@inline Base.@propagate_inbounds function heat_cap_h2_kernel(mc::DQMC, model::TwoBandModel, sites::NTuple{4, Int}, 
    G::GreensMatrix, flv) 
    return heat_cap_h2_kernel(mc, model, sites, G, flv, field(mc))
end


@inline Base.@propagate_inbounds function heat_cap_h2_kernel(
        mc, model::TwoBandModel, sites::NTuple{4, Int}, packed_greens::GreensMatrix, 
        flv, ::AbstractMagnBosonField  )
    i, iPb1, j, jPb2 = sites  #i, i+b₁, j, j+b₂
	G = packed_greens
    N = length(lattice(model))
    T = mc.stack.hopping_matrix

    return (G.val[i, iPb1]*T[iPb1, i] + G.val[i + N, iPb1 + N]*T[iPb1 + N, i + N] + G.val[i + 2N, iPb1 + 2N]*T[iPb1 + 2N, i + 2N] + 
    G.val[i + 3N, iPb1 + 3N]*T[iPb1 + 3N, i + 3N] - I[i, iPb1]*(T[iPb1, i] + T[iPb1 + N, i + N] + T[iPb1 + 2N, i + 2N] + T[iPb1 + 3N, i + 3N]))*
   (G.val[j, jPb2]*T[jPb2, j] + G.val[j + N, jPb2 + N]*T[jPb2 + N, j + N] + G.val[j + 2N, jPb2 + 2N]*T[jPb2 + 2N, j + 2N] + 
    G.val[j + 3N, jPb2 + 3N]*T[jPb2 + 3N, j + 3N] - I[j, jPb2]*(T[jPb2, j] + T[jPb2 + N, j + N] + T[jPb2 + 2N, j + 2N] + T[jPb2 + 3N, j + 3N]))+
    G.val[i, jPb2]*(-G.val[j, iPb1] + I[j, iPb1])*T[iPb1, i]*T[jPb2, j] - G.val[j, iPb1 + N]*G.val[i + N, jPb2]*T[jPb2, j]*T[iPb1 + N, i + N] - 
   G.val[i, jPb2 + N]*G.val[j + N, iPb1]*T[iPb1, i]*T[jPb2 + N, j + N] + G.val[i + N, jPb2 + N]*(-G.val[j + N, iPb1 + N] + I[j, iPb1])*T[iPb1 + N, i + N]*
    T[jPb2 + N, j + N] - G.val[j, iPb1 + 2N]*G.val[i + 2N, jPb2]*T[jPb2, j]*T[iPb1 + 2N, i + 2N] - G.val[j + N, iPb1 + 2N]*G.val[i + 2N, jPb2 + N]*T[jPb2 + N, j + N]*
    T[iPb1 + 2N, i + 2N] - G.val[i, jPb2 + 2N]*G.val[j + 2N, iPb1]*T[iPb1, i]*T[jPb2 + 2N, j + 2N] - 
   G.val[i + N, jPb2 + 2N]*G.val[j + 2N, iPb1 + N]*T[iPb1 + N, i + N]*T[jPb2 + 2N, j + 2N] + G.val[i + 2N, jPb2 + 2N]*(-G.val[j + 2N, iPb1 + 2N] + I[j, iPb1])*
    T[iPb1 + 2N, i + 2N]*T[jPb2 + 2N, j + 2N] - G.val[j, iPb1 + 3N]*G.val[i + 3N, jPb2]*T[jPb2, j]*T[iPb1 + 3N, i + 3N] - 
   G.val[j + N, iPb1 + 3N]*G.val[i + 3N, jPb2 + N]*T[jPb2 + N, j + N]*T[iPb1 + 3N, i + 3N] - G.val[j + 2N, iPb1 + 3N]*G.val[i + 3N, jPb2 + 2N]*T[jPb2 + 2N, j + 2N]*
    T[iPb1 + 3N, i + 3N] - G.val[i, jPb2 + 3N]*G.val[j + 3N, iPb1]*T[iPb1, i]*T[jPb2 + 3N, j + 3N] - 
   G.val[i + N, jPb2 + 3N]*G.val[j + 3N, iPb1 + N]*T[iPb1 + N, i + N]*T[jPb2 + 3N, j + 3N] - G.val[i + 2N, jPb2 + 3N]*G.val[j + 3N, iPb1 + 2N]*T[iPb1 + 2N, i + 2N]*
    T[jPb2 + 3N, j + 3N] + G.val[i + 3N, jPb2 + 3N]*(-G.val[j + 3N, iPb1 + 3N] + I[j, iPb1])*T[iPb1 + 3N, i + 3N]*T[jPb2 + 3N, j + 3N]
end
@inline Base.@propagate_inbounds function heat_cap_h2_kernel(
    mc, model::TwoBandModel, sites::NTuple{4, Int}, packed_greens::GreensMatrix, 
    flv, ::Union{Discrete_MBF1_X_symm, Discrete_MBF2_symm} )
    i, iPb1, j, jPb2 = sites  #i, i+b₁, j, j+b₂
	G = packed_greens
    N = length(lattice(model))
    T = mc.stack.hopping_matrix

    return 4*real(G.val[i, iPb1]*T[iPb1, i] + G.val[i + N, iPb1 + N]*T[iPb1 + N, i + N] - I[i, iPb1]*(T[iPb1, i] + T[iPb1 + N, i + N]))*
    real(G.val[j, jPb2]*T[jPb2, j] + G.val[j + N, jPb2 + N]*T[jPb2 + N, j + N] - I[j, jPb2]*(T[jPb2, j] + T[jPb2 + N, j + N]))+
    2*real(G.val[i, jPb2]*(-G.val[j, iPb1] + I[j, iPb1])*T[iPb1, i]*T[jPb2, j] - G.val[j, iPb1 + N]*G.val[i + N, jPb2]*T[jPb2, j]*T[iPb1 + N, i + N] - 
      (G.val[i, jPb2 + N]*G.val[j + N, iPb1]*T[iPb1, i] + G.val[i + N, jPb2 + N]*(G.val[j + N, iPb1 + N] - I[j, iPb1])*T[iPb1 + N, i + N])*T[jPb2 + N, j + N])

end


###########################
### heat capacity kernel for h₃
###########################


"""
Calculates the heat capacity kernel for h₃
"""
@inline Base.@propagate_inbounds function heat_cap_h3_kernel(mc::DQMC, model::TwoBandModel, sites::NTuple{3, Int}, 
    G::GreensMatrix, flv) 
    return heat_cap_h3_kernel(mc, model, sites, G, flv, field(mc))
end




@inline Base.@propagate_inbounds function heat_cap_h3_kernel(
        mc, model::TwoBandModel, sites::NTuple{3, Int}, packed_greens::GreensMatrix, 
        flv, ::Discrete_MBF1_X  )
    i, iPb, j = sites  #i, i+b, j
    G = packed_greens
    N = length(lattice(model))
    T = mc.stack.hopping_matrix

    h3a=-2*(-2*G.val[j, j + 2N]*G.val[j + N, j + 3N] + G.val[j + N, j + N]*(1 - 2*G.val[j + 2N, j + 2N]) + G.val[j + 2N, j + 2N] + 
    2*(G.val[j + 2N, j + N]*G.val[j + 3N, j] + G.val[j + N, j + 2N]*(G.val[j + 2N, j + N] + G.val[j + 3N, j]) + 
      G.val[j, j + 3N]*(G.val[j + N, j + 2N] + G.val[j + 2N, j + N] + G.val[j + 3N, j])) - 
    2*(G.val[j, j + N]*G.val[j + 2N, j + 3N] + G.val[j + 2N, j]*G.val[j + 3N, j + N] + G.val[j + N, j]*G.val[j + 3N, j + 2N]) + G.val[j, j]*(1 - 2*G.val[j + 3N, j + 3N]) + 
    G.val[j + 3N, j + 3N])*(G.val[i, iPb]*T[iPb, i] + G.val[i + N, iPb + N]*T[iPb + N, i + N] + G.val[i + 2N, iPb + 2N]*T[iPb + 2N, i + 2N] + 
    G.val[i + 3N, iPb + 3N]*T[iPb + 3N, i + 3N] - I[i, iPb]*(T[iPb, i] + T[iPb + N, i + N] + T[iPb + 2N, i + 2N] + T[iPb + 3N, i + 3N]));

    h3b=2*(G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 2N, j + N] + G.val[j + 3N, j])*(2*G.val[i, j + 2N]*G.val[j + N, iPb]*T[iPb, i] + 
    2*G.val[i, j + N]*G.val[j + 2N, iPb]*T[iPb, i] + 2*G.val[i, j]*G.val[j + 3N, iPb]*T[iPb, i] - G.val[j + 3N, iPb]*I[i, j]*T[iPb, i] + 
    G.val[i, j + 3N]*(2*G.val[j, iPb] - I[j, iPb])*T[iPb, i] - G.val[j + 2N, iPb + N]*I[i, j]*T[iPb + N, i + N] - G.val[i + N, j + 2N]*I[j, iPb]*T[iPb + N, i + N] - 
    G.val[j + N, iPb + 2N]*I[i, j]*T[iPb + 2N, i + 2N] - G.val[i + 2N, j + N]*I[j, iPb]*T[iPb + 2N, i + 2N] - 
    (G.val[j, iPb + 3N]*I[i, j] + G.val[i + 3N, j]*I[j, iPb])*T[iPb + 3N, i + 3N] + 
    2*((G.val[j, iPb + N]*G.val[i + N, j + 3N] + G.val[i + N, j + 2N]*G.val[j + N, iPb + N] + G.val[i + N, j + N]*G.val[j + 2N, iPb + N] + G.val[i + N, j]*G.val[j + 3N, iPb + N])*
       T[iPb + N, i + N] + (G.val[j + N, iPb + 2N]*G.val[i + 2N, j + 2N] + G.val[j, iPb + 2N]*G.val[i + 2N, j + 3N] + G.val[i + 2N, j + N]*G.val[j + 2N, iPb + 2N] + 
        G.val[i + 2N, j]*G.val[j + 3N, iPb + 2N])*T[iPb + 2N, i + 2N] + (G.val[j + 2N, iPb + 3N]*G.val[i + 3N, j + N] + G.val[j + N, iPb + 3N]*G.val[i + 3N, j + 2N] + 
        G.val[j, iPb + 3N]*G.val[i + 3N, j + 3N] + G.val[i + 3N, j]*G.val[j + 3N, iPb + 3N])*T[iPb + 3N, i + 3N]));

    h3c=-4*G.val[i, j + 3N]*G.val[j, j + 2N]*G.val[j + N, iPb]*T[iPb, i] - 4*G.val[i, j + 2N]*G.val[j + N, iPb]*G.val[j + N, j + 2N]*T[iPb, i] - 
    4*G.val[i, j + 3N]*G.val[j, j + N]*G.val[j + 2N, iPb]*T[iPb, i] + 2*G.val[i, j + 2N]*(1 - 2*G.val[j + N, j + N])*G.val[j + 2N, iPb]*T[iPb, i] - 
    4*G.val[i, j + N]*G.val[j + 2N, iPb]*G.val[j + 2N, j + N]*T[iPb, i] + 2*G.val[i, j + N]*G.val[j + N, iPb]*(1 - 2*G.val[j + 2N, j + 2N])*T[iPb, i] + 
    2*G.val[i, j + 3N]*(1 - 2*G.val[j, j])*G.val[j + 3N, iPb]*T[iPb, i] - 4*G.val[i, j + 2N]*G.val[j + N, j]*G.val[j + 3N, iPb]*T[iPb, i] - 
    4*G.val[i, j + N]*G.val[j + 2N, j]*G.val[j + 3N, iPb]*T[iPb, i] + 2*G.val[j + 3N, iPb]*G.val[j + 3N, j]*(-2*G.val[i, j] + I[i, j])*T[iPb, i] + 
    2*G.val[j + 2N, iPb]*G.val[j + 3N, j + N]*(-2*G.val[i, j] + I[i, j])*T[iPb, i] + 2*G.val[j + N, iPb]*G.val[j + 3N, j + 2N]*(-2*G.val[i, j] + I[i, j])*T[iPb, i] + 
    2*G.val[i, j + 3N]*G.val[j, j + 3N]*(-2*G.val[j, iPb] + I[j, iPb])*T[iPb, i] + 2*G.val[i, j + 2N]*G.val[j + N, j + 3N]*(-2*G.val[j, iPb] + I[j, iPb])*T[iPb, i] + 
    2*G.val[i, j + N]*G.val[j + 2N, j + 3N]*(-2*G.val[j, iPb] + I[j, iPb])*T[iPb, i] + 
    (-1 + 2*G.val[j + 3N, j + 3N])*(G.val[j, iPb]*I[i, j] + G.val[i, j]*(-2*G.val[j, iPb] + I[j, iPb]))*T[iPb, i] - 
    4*G.val[j, iPb + N]*G.val[j, j + 3N]*G.val[i + N, j + 3N]*T[iPb + N, i + N] - 4*G.val[j, iPb + N]*G.val[i + N, j + 2N]*G.val[j + N, j + 3N]*T[iPb + N, i + N] - 
    4*G.val[j, j + N]*G.val[i + N, j + 3N]*G.val[j + 2N, iPb + N]*T[iPb + N, i + N] + 2*G.val[i + N, j + 2N]*(1 - 2*G.val[j + N, j + N])*G.val[j + 2N, iPb + N]*
     T[iPb + N, i + N] + 2*(1 - 2*G.val[j, j])*G.val[i + N, j + 3N]*G.val[j + 3N, iPb + N]*T[iPb + N, i + N] - 
    4*G.val[i + N, j + 2N]*G.val[j + N, j]*G.val[j + 3N, iPb + N]*T[iPb + N, i + N] - 4*G.val[i + N, j]*G.val[j + 3N, j]*G.val[j + 3N, iPb + N]*T[iPb + N, i + N] - 
    4*G.val[i + N, j]*G.val[j + 2N, iPb + N]*G.val[j + 3N, j + N]*T[iPb + N, i + N] + 2*G.val[j, iPb + N]*G.val[i + N, j]*(1 - 2*G.val[j + 3N, j + 3N])*T[iPb + N, i + N] + 
    2*G.val[j + 2N, iPb + N]*G.val[j + 2N, j + N]*(-2*G.val[i + N, j + N] + I[i, j])*T[iPb + N, i + N] + 
    2*G.val[j, iPb + N]*G.val[j + 2N, j + 3N]*(-2*G.val[i + N, j + N] + I[i, j])*T[iPb + N, i + N] + 
    2*G.val[j + 2N, j]*G.val[j + 3N, iPb + N]*(-2*G.val[i + N, j + N] + I[i, j])*T[iPb + N, i + N] + 
    2*G.val[j, j + 2N]*G.val[i + N, j + 3N]*(-2*G.val[j + N, iPb + N] + I[j, iPb])*T[iPb + N, i + N] + 
    2*G.val[i + N, j + 2N]*G.val[j + N, j + 2N]*(-2*G.val[j + N, iPb + N] + I[j, iPb])*T[iPb + N, i + N] + 
    2*G.val[i + N, j]*G.val[j + 3N, j + 2N]*(-2*G.val[j + N, iPb + N] + I[j, iPb])*T[iPb + N, i + N] + 
    (-1 + 2*G.val[j + 2N, j + 2N])*(G.val[j + N, iPb + N]*I[i, j] + G.val[i + N, j + N]*(-2*G.val[j + N, iPb + N] + I[j, iPb]))*T[iPb + N, i + N] - 
    4*G.val[j, iPb + 2N]*G.val[j, j + 3N]*G.val[i + 2N, j + 3N]*T[iPb + 2N, i + 2N] - 4*G.val[j, j + 2N]*G.val[j + N, iPb + 2N]*G.val[i + 2N, j + 3N]*
     T[iPb + 2N, i + 2N] + 2*G.val[j + N, iPb + 2N]*G.val[i + 2N, j + N]*(1 - 2*G.val[j + 2N, j + 2N])*T[iPb + 2N, i + 2N] - 
    4*G.val[j, iPb + 2N]*G.val[i + 2N, j + N]*G.val[j + 2N, j + 3N]*T[iPb + 2N, i + 2N] + 2*(1 - 2*G.val[j, j])*G.val[i + 2N, j + 3N]*G.val[j + 3N, iPb + 2N]*
     T[iPb + 2N, i + 2N] - 4*G.val[i + 2N, j + N]*G.val[j + 2N, j]*G.val[j + 3N, iPb + 2N]*T[iPb + 2N, i + 2N] - 
    4*G.val[i + 2N, j]*G.val[j + 3N, j]*G.val[j + 3N, iPb + 2N]*T[iPb + 2N, i + 2N] - 4*G.val[j + N, iPb + 2N]*G.val[i + 2N, j]*G.val[j + 3N, j + 2N]*
     T[iPb + 2N, i + 2N] + 2*G.val[j, iPb + 2N]*G.val[i + 2N, j]*(1 - 2*G.val[j + 3N, j + 3N])*T[iPb + 2N, i + 2N] + 
    2*G.val[j + N, iPb + 2N]*G.val[j + N, j + 2N]*(-2*G.val[i + 2N, j + 2N] + I[i, j])*T[iPb + 2N, i + 2N] + 
    2*G.val[j, iPb + 2N]*G.val[j + N, j + 3N]*(-2*G.val[i + 2N, j + 2N] + I[i, j])*T[iPb + 2N, i + 2N] + 
    2*G.val[j + N, j]*G.val[j + 3N, iPb + 2N]*(-2*G.val[i + 2N, j + 2N] + I[i, j])*T[iPb + 2N, i + 2N] + 
    2*G.val[j, j + N]*G.val[i + 2N, j + 3N]*(-2*G.val[j + 2N, iPb + 2N] + I[j, iPb])*T[iPb + 2N, i + 2N] + 
    2*G.val[i + 2N, j + N]*G.val[j + 2N, j + N]*(-2*G.val[j + 2N, iPb + 2N] + I[j, iPb])*T[iPb + 2N, i + 2N] + 
    2*G.val[i + 2N, j]*G.val[j + 3N, j + N]*(-2*G.val[j + 2N, iPb + 2N] + I[j, iPb])*T[iPb + 2N, i + 2N] + 
    (-1 + 2*G.val[j + N, j + N])*(G.val[j + 2N, iPb + 2N]*I[i, j] + G.val[i + 2N, j + 2N]*(-2*G.val[j + 2N, iPb + 2N] + I[j, iPb]))*T[iPb + 2N, i + 2N] - 
    (4*G.val[j + 2N, j + N]*G.val[j + 2N, iPb + 3N]*G.val[i + 3N, j + N] + G.val[j + N, iPb + 3N]*((-2 + 4*G.val[j + 2N, j + 2N])*G.val[i + 3N, j + N] + 
        4*G.val[j + N, j + 2N]*G.val[i + 3N, j + 2N]) + 4*G.val[j, j + 2N]*G.val[j + N, iPb + 3N]*G.val[i + 3N, j + 3N] + 
      4*G.val[j + N, iPb + 3N]*G.val[i + 3N, j]*G.val[j + 3N, j + 2N] - 2*G.val[i + 3N, j + 3N]*G.val[j + 3N, iPb + 3N] + 
      4*G.val[j, j]*G.val[i + 3N, j + 3N]*G.val[j + 3N, iPb + 3N] + 4*G.val[i + 3N, j]*G.val[j + 3N, j]*G.val[j + 3N, iPb + 3N] + 
      G.val[j + 2N, iPb + 3N]*((-2 + 4*G.val[j + N, j + N])*G.val[i + 3N, j + 2N] + 4*G.val[i + 3N, j]*G.val[j + 3N, j + N] + 
        G.val[j, j + N]*(4*G.val[i + 3N, j + 3N] - 2*I[i, j])) + G.val[j, iPb + 3N]*(4*G.val[j + 2N, j + 3N]*G.val[i + 3N, j + N] + 
        4*G.val[j + N, j + 3N]*G.val[i + 3N, j + 2N] + G.val[i + 3N, j]*(-2 + 4*G.val[j + 3N, j + 3N]) + 2*G.val[j, j + 3N]*(2*G.val[i + 3N, j + 3N] - I[i, j])) - 
      2*G.val[j, j + 2N]*G.val[j + N, iPb + 3N]*I[i, j] + G.val[j + 3N, iPb + 3N]*I[i, j] - 2*G.val[j, j]*G.val[j + 3N, iPb + 3N]*I[i, j] + 
      2*(G.val[j + 2N, j]*G.val[i + 3N, j + N] + G.val[j + N, j]*G.val[i + 3N, j + 2N])*(2*G.val[j + 3N, iPb + 3N] - I[j, iPb]) + 
      ((1 - 2*G.val[j, j])*G.val[i + 3N, j + 3N] - 2*G.val[i + 3N, j]*G.val[j + 3N, j])*I[j, iPb])*T[iPb + 3N, i + 3N];

    
    return h3a + h3b +h3c 
end
@inline Base.@propagate_inbounds function heat_cap_h3_kernel(
    mc, model::TwoBandModel, sites::NTuple{3, Int}, packed_greens::GreensMatrix, 
    flv, ::Discrete_MBF1_X_symm  )
    i, iPb, j = sites  #i, i+b, j
    G = packed_greens
    N = length(lattice(model))
    T = mc.stack.hopping_matrix

    h3a=-8*(abs2(G.val[j, j + N]) + abs2(G.val[j + N, j]) + 
      real(G.val[j, j] + 2*conj(G.val[j + N, j])*G.val[j, j + N] + 2*G.val[j, j + N]*G.val[j + N, j] + 
      G.val[j + N, j + N] - 2*G.val[j, j]*G.val[j + N, j + N]))*
      real(G.val[i, iPb]*T[iPb, i] + G.val[i + N, iPb + N]*T[iPb + N, i + N] - I[i, iPb]*(T[iPb, i] + T[iPb + N, i + N]));

    h3b=8*real(G.val[j, j + N] + G.val[j + N, j])*real(2*(G.val[i, j + N]*G.val[j, iPb] + G.val[i, j]*G.val[j + N, iPb])*T[iPb, i] - G.val[j + N, iPb]*I[i, j]*T[iPb, i] - 
    G.val[i, j + N]*I[j, iPb]*T[iPb, i] + 2*(G.val[j, iPb + N]*G.val[i + N, j + N] + G.val[i + N, j]*G.val[j + N, iPb + N])*T[iPb + N, i + N] - 
    G.val[j, iPb + N]*I[i, j]*T[iPb + N, i + N] - G.val[i + N, j]*I[j, iPb]*T[iPb + N, i + N]);

    h3c=2*real(2*G.val[i, j + N]*(1 - 2*G.val[j, j])*G.val[j + N, iPb]*T[iPb, i] + 2*G.val[j + N, iPb]*G.val[j + N, j]*(-2*G.val[i, j] + I[i, j])*T[iPb, i] + 
    2*G.val[i, j + N]*G.val[j, j + N]*(-2*G.val[j, iPb] + I[j, iPb])*T[iPb, i] + (-1 + 2*G.val[j + N, j + N])*
     (G.val[j, iPb]*I[i, j] + G.val[i, j]*(-2*G.val[j, iPb] + I[j, iPb]))*T[iPb, i] + 2*G.val[j, iPb + N]*G.val[i + N, j]*(1 - 2*G.val[j + N, j + N])*T[iPb + N, i + N] + 
    2*G.val[j, iPb + N]*G.val[j, j + N]*(-2*G.val[i + N, j + N] + I[i, j])*T[iPb + N, i + N] + 2*G.val[i + N, j]*G.val[j + N, j]*(-2*G.val[j + N, iPb + N] + I[j, iPb])*
     T[iPb + N, i + N] + (-1 + 2*G.val[j, j])*(G.val[j + N, iPb + N]*I[i, j] + G.val[i + N, j + N]*(-2*G.val[j + N, iPb + N] + I[j, iPb]))*T[iPb + N, i + N]);

    return h3a + h3b +h3c 
end

###########################
### heat capacity kernel for h₄
###########################


"""
Calculates the heat capacity kernel for h₄
"""
@inline Base.@propagate_inbounds function heat_cap_h4_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv) 
    return heat_cap_h4_kernel(mc, model, ij, G, flv, field(mc))
end




@inline Base.@propagate_inbounds function heat_cap_h4_kernel(
        mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::GreensMatrix, 
        flv, ::Discrete_MBF1_X  )
    i, j = ij   
	G = packed_greens
    N = length(lattice(model))

    h4a=(-2*G.val[i, i + 2N]*G.val[i + N, i + 3N] + G.val[i + N, i + N]*(1 - 2*G.val[i + 2N, i + 2N]) + G.val[i + 2N, i + 2N] + 
    2*(G.val[i + 2N, i + N]*G.val[i + 3N, i] + G.val[i + N, i + 2N]*(G.val[i + 2N, i + N] + G.val[i + 3N, i]) + 
      G.val[i, i + 3N]*(G.val[i + N, i + 2N] + G.val[i + 2N, i + N] + G.val[i + 3N, i])) - 
    2*(G.val[i, i + N]*G.val[i + 2N, i + 3N] + G.val[i + 2N, i]*G.val[i + 3N, i + N] + G.val[i + N, i]*G.val[i + 3N, i + 2N]) + G.val[i, i]*(1 - 2*G.val[i + 3N, i + 3N]) + 
    G.val[i + 3N, i + 3N])*(-2*G.val[j, j + 2N]*G.val[j + N, j + 3N] + G.val[j + N, j + N]*(1 - 2*G.val[j + 2N, j + 2N]) + G.val[j + 2N, j + 2N] + 
    2*(G.val[j + 2N, j + N]*G.val[j + 3N, j] + G.val[j + N, j + 2N]*(G.val[j + 2N, j + N] + G.val[j + 3N, j]) + 
      G.val[j, j + 3N]*(G.val[j + N, j + 2N] + G.val[j + 2N, j + N] + G.val[j + 3N, j])) - 
    2*(G.val[j, j + N]*G.val[j + 2N, j + 3N] + G.val[j + 2N, j]*G.val[j + 3N, j + N] + G.val[j + N, j]*G.val[j + 3N, j + 2N]) + G.val[j, j]*(1 - 2*G.val[j + 3N, j + 3N]) + 
    G.val[j + 3N, j + 3N]);

    h4b=4*((G.val[i, j + 3N]*G.val[i + N, j + 2N] - G.val[i, j + 2N]*G.val[i + N, j + 3N])*(G.val[j, i + 3N]*G.val[j + N, i + 2N] - G.val[j, i + 2N]*G.val[j + N, i + 3N]) + 
    (G.val[i, j + 3N]*G.val[i + 2N, j + N] - G.val[i, j + N]*G.val[i + 2N, j + 3N])*(G.val[j, i + 3N]*G.val[j + 2N, i + N] - G.val[j, i + N]*G.val[j + 2N, i + 3N]) + 
    (G.val[j + N, i + 3N]*G.val[j + 2N, i] - G.val[j + N, i]*G.val[j + 2N, i + 3N])*(G.val[i, j + 2N]*G.val[i + 3N, j + N] - G.val[i, j + N]*G.val[i + 3N, j + 2N]) + 
    (G.val[i + 2N, j + N]*G.val[i + 3N, j] - G.val[i + 2N, j]*G.val[i + 3N, j + N])*(G.val[j + 2N, i + N]*G.val[j + 3N, i] - G.val[j + 2N, i]*G.val[j + 3N, i + N]) + 
    (G.val[i + N, j + 3N]*G.val[i + 2N, j] - G.val[i + N, j]*G.val[i + 2N, j + 3N])*(G.val[j, i + 2N]*G.val[j + 3N, i + N] - G.val[j, i + N]*G.val[j + 3N, i + 2N]) + 
    (G.val[i + N, j + 2N]*G.val[i + 3N, j] - G.val[i + N, j]*G.val[i + 3N, j + 2N])*(G.val[j + N, i + 2N]*G.val[j + 3N, i] - G.val[j + N, i]*G.val[j + 3N, i + 2N]) + 
    (G.val[i + N, j + 3N]*G.val[i + 2N, j + 2N] - G.val[i + N, j + 2N]*G.val[i + 2N, j + 3N])*(-(G.val[j, i + N]*G.val[j + N, i + 2N]) + 
      G.val[j, i + 2N]*(G.val[j + N, i + N] - I[j, i])) + (G.val[i, j + 3N]*G.val[i + 2N, j + 2N] - G.val[i, j + 2N]*G.val[i + 2N, j + 3N])*
     (-(G.val[j, i + N]*G.val[j + N, i + 3N]) + G.val[j, i + 3N]*(G.val[j + N, i + N] - I[j, i])) + 
    (G.val[i + 2N, j + 2N]*G.val[i + 3N, j + N] - G.val[i + 2N, j + N]*G.val[i + 3N, j + 2N])*(-(G.val[j + N, i]*G.val[j + 2N, i + N]) + 
      G.val[j + 2N, i]*(G.val[j + N, i + N] - I[j, i])) + (G.val[i + 2N, j + 2N]*G.val[i + 3N, j] - G.val[i + 2N, j]*G.val[i + 3N, j + 2N])*
     (-(G.val[j + N, i]*G.val[j + 3N, i + N]) + G.val[j + 3N, i]*(G.val[j + N, i + N] - I[j, i])) + 
    (G.val[i, j + 3N]*G.val[i + N, j + N] - G.val[i, j + N]*G.val[i + N, j + 3N])*(-(G.val[j, i + 2N]*G.val[j + 2N, i + 3N]) + 
      G.val[j, i + 3N]*(G.val[j + 2N, i + 2N] - I[j, i])) + (G.val[i, j + 2N]*G.val[i + N, j + N] - G.val[i, j + N]*G.val[i + N, j + 2N])*
     (-(G.val[j + N, i + 2N]*G.val[j + 2N, i + 3N]) + G.val[j + N, i + 3N]*(G.val[j + 2N, i + 2N] - I[j, i])) + 
    (G.val[i + N, j + N]*G.val[i + 3N, j] - G.val[i + N, j]*G.val[i + 3N, j + N])*(-(G.val[j + 2N, i]*G.val[j + 3N, i + 2N]) + 
      G.val[j + 3N, i]*(G.val[j + 2N, i + 2N] - I[j, i])) + (G.val[i + N, j + N]*G.val[i + 2N, j] - G.val[i + N, j]*G.val[i + 2N, j + N])*
     (-(G.val[j + 2N, i + N]*G.val[j + 3N, i + 2N]) + G.val[j + 3N, i + N]*(G.val[j + 2N, i + 2N] - I[j, i])) + 
    (G.val[i + N, j + 3N]*G.val[i + 3N, j + 2N] - G.val[i + N, j + 2N]*G.val[i + 3N, j + 3N])*(G.val[j, i + 2N]*G.val[j + N, i] + 
      G.val[j + N, i + 2N]*(-G.val[j, i] + I[j, i])) + (G.val[i, j + 3N]*G.val[i + 3N, j + 2N] - G.val[i, j + 2N]*G.val[i + 3N, j + 3N])*
     (G.val[j, i + 3N]*G.val[j + N, i] + G.val[j + N, i + 3N]*(-G.val[j, i] + I[j, i])) + 
    (G.val[i + 2N, j + 3N]*G.val[i + 3N, j + N] - G.val[i + 2N, j + N]*G.val[i + 3N, j + 3N])*(G.val[j, i + N]*G.val[j + 2N, i] + 
      G.val[j + 2N, i + N]*(-G.val[j, i] + I[j, i])) + (G.val[i, j + 3N]*G.val[i + 3N, j + N] - G.val[i, j + N]*G.val[i + 3N, j + 3N])*
     (G.val[j, i + 3N]*G.val[j + 2N, i] + G.val[j + 2N, i + 3N]*(-G.val[j, i] + I[j, i])) + 
    (G.val[i + 2N, j + 3N]*G.val[i + 3N, j] - G.val[i + 2N, j]*G.val[i + 3N, j + 3N])*(G.val[j, i + N]*G.val[j + 3N, i] + G.val[j + 3N, i + N]*(-G.val[j, i] + I[j, i])) + 
    (G.val[i + N, j + 3N]*G.val[i + 3N, j] - G.val[i + N, j]*G.val[i + 3N, j + 3N])*(G.val[j, i + 2N]*G.val[j + 3N, i] + G.val[j + 3N, i + 2N]*(-G.val[j, i] + I[j, i])) + 
    (G.val[i, j + 2N]*G.val[i + 2N, j + N] - G.val[i, j + N]*G.val[i + 2N, j + 2N])*(G.val[j + N, i + 3N]*G.val[j + 2N, i + N] + 
      G.val[j + 2N, i + 3N]*(-G.val[j + N, i + N] + I[j, i])) + (G.val[i + N, j + 2N]*G.val[i + 2N, j] - G.val[i + N, j]*G.val[i + 2N, j + 2N])*
     (G.val[j + N, i + 2N]*G.val[j + 3N, i + N] + G.val[j + 3N, i + 2N]*(-G.val[j + N, i + N] + I[j, i])) + 
    (G.val[i + 2N, j + 3N]*G.val[i + 3N, j + 2N] - G.val[i + 2N, j + 2N]*G.val[i + 3N, j + 3N])*
     (G.val[j, i + N]*G.val[j + N, i] + (G.val[j, i] - I[j, i])*(-G.val[j + N, i + N] + I[j, i])) + 
    (G.val[i + N, j + 3N]*G.val[i + 2N, j + N] - G.val[i + N, j + N]*G.val[i + 2N, j + 3N])*(G.val[j, i + 2N]*G.val[j + 2N, i + N] + 
      G.val[j, i + N]*(-G.val[j + 2N, i + 2N] + I[j, i])) + (G.val[i + N, j + 2N]*G.val[i + 3N, j + N] - G.val[i + N, j + N]*G.val[i + 3N, j + 2N])*
     (G.val[j + N, i + 2N]*G.val[j + 2N, i] + G.val[j + N, i]*(-G.val[j + 2N, i + 2N] + I[j, i])) + 
    (G.val[i + N, j + 3N]*G.val[i + 3N, j + N] - G.val[i + N, j + N]*G.val[i + 3N, j + 3N])*(G.val[j, i + 2N]*G.val[j + 2N, i] + 
      (G.val[j, i] - I[j, i])*(-G.val[j + 2N, i + 2N] + I[j, i])) + (G.val[i + N, j + 2N]*G.val[i + 2N, j + N] - G.val[i + N, j + N]*G.val[i + 2N, j + 2N])*
     (G.val[j + N, i + 2N]*G.val[j + 2N, i + N] + (G.val[j + N, i + N] - I[j, i])*(-G.val[j + 2N, i + 2N] + I[j, i])) + 
    (G.val[i, j + 3N]*G.val[i + 2N, j] - G.val[i, j]*G.val[i + 2N, j + 3N])*(G.val[j, i + 3N]*G.val[j + 3N, i + N] + G.val[j, i + N]*(-G.val[j + 3N, i + 3N] + I[j, i])) + 
    (G.val[i, j + 3N]*G.val[i + N, j] - G.val[i, j]*G.val[i + N, j + 3N])*(G.val[j, i + 3N]*G.val[j + 3N, i + 2N] + G.val[j, i + 2N]*(-G.val[j + 3N, i + 3N] + I[j, i])) + 
    (G.val[i, j + 2N]*G.val[i + 3N, j] - G.val[i, j]*G.val[i + 3N, j + 2N])*(G.val[j + N, i + 3N]*G.val[j + 3N, i] + G.val[j + N, i]*(-G.val[j + 3N, i + 3N] + I[j, i])) + 
    (G.val[i, j + 2N]*G.val[i + N, j] - G.val[i, j]*G.val[i + N, j + 2N])*(G.val[j + N, i + 3N]*G.val[j + 3N, i + 2N] + 
      G.val[j + N, i + 2N]*(-G.val[j + 3N, i + 3N] + I[j, i])) + (G.val[i, j + N]*G.val[i + 3N, j] - G.val[i, j]*G.val[i + 3N, j + N])*
     (G.val[j + 2N, i + 3N]*G.val[j + 3N, i] + G.val[j + 2N, i]*(-G.val[j + 3N, i + 3N] + I[j, i])) + (G.val[i, j + N]*G.val[i + 2N, j] - G.val[i, j]*G.val[i + 2N, j + N])*
     (G.val[j + 2N, i + 3N]*G.val[j + 3N, i + N] + G.val[j + 2N, i + N]*(-G.val[j + 3N, i + 3N] + I[j, i])) + 
    (G.val[i, j + 3N]*G.val[i + 3N, j] - G.val[i, j]*G.val[i + 3N, j + 3N])*(G.val[j, i + 3N]*G.val[j + 3N, i] + 
      (G.val[j, i] - I[j, i])*(-G.val[j + 3N, i + 3N] + I[j, i])) + (G.val[i, j + 2N]*G.val[i + 2N, j] - G.val[i, j]*G.val[i + 2N, j + 2N])*
     (G.val[j + N, i + 3N]*G.val[j + 3N, i + N] + (G.val[j + N, i + N] - I[j, i])*(-G.val[j + 3N, i + 3N] + I[j, i])) + 
    (G.val[i, j + N]*G.val[i + N, j] - G.val[i, j]*G.val[i + N, j + N])*(G.val[j + 2N, i + 3N]*G.val[j + 3N, i + 2N] + 
      (G.val[j + 2N, i + 2N] - I[j, i])*(-G.val[j + 3N, i + 3N] + I[j, i])));

    h4c=(4*(G.val[i + N, j + 2N]*(-1/2 + G.val[i + 2N, i + 2N]) + G.val[i, j + 2N]*G.val[i + 2N, i + 3N] - 
    G.val[i + 2N, j + 2N]*(G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 3N, i]) + G.val[i + 2N, i]*G.val[i + 3N, j + 2N])*
   (-(G.val[j, i + N]*G.val[j + N, j + 3N]) + (1/2 - G.val[j + N, j + N])*G.val[j + 2N, i + N] - G.val[j + N, j]*G.val[j + 3N, i + N] + 
    (G.val[j, j + 3N] + G.val[j + 2N, j + N] + G.val[j + 3N, j])*(G.val[j + N, i + N] - I[j, i])) - 
  (2*G.val[i, j + 3N]*G.val[i + N, i + 3N] + (-1 + 2*G.val[i + N, i + N])*G.val[i + 2N, j + 3N] - 
    2*G.val[i + N, j + 3N]*(G.val[i, i + 3N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) + 2*G.val[i + N, i]*G.val[i + 3N, j + 3N])*
   (2*G.val[j, j + 2N]*G.val[j + N, i + 2N] - 2*G.val[j, i + 2N]*(G.val[j + N, j + 2N] + G.val[j + 2N, j + N] + G.val[j + 3N, j]) + 
    (-1 + 2*G.val[j, j])*G.val[j + 3N, i + 2N] + 2*G.val[j, j + N]*(G.val[j + 2N, i + 2N] - I[j, i])) + 
  (-2*G.val[i, j + N]*G.val[i + N, i + 3N] + (1 - 2*G.val[i + N, i + N])*G.val[i + 2N, j + N] + 
    2*G.val[i + N, j + N]*(G.val[i, i + 3N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) - 2*G.val[i + N, i]*G.val[i + 3N, j + N])*
   (G.val[j + N, i + 2N]*(-1 + 2*G.val[j + 2N, j + 2N]) + 2*G.val[j, i + 2N]*G.val[j + 2N, j + 3N] + 2*G.val[j + 2N, j]*G.val[j + 3N, i + 2N] - 
    2*(G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 3N, j])*(G.val[j + 2N, i + 2N] - I[j, i])) - 
  (-2*G.val[i, i + 2N]*G.val[i + N, j + 3N] - 2*G.val[i, i + N]*G.val[i + 2N, j + 3N] + 
    2*G.val[i, j + 3N]*(G.val[i + N, i + 2N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) + (1 - 2*G.val[i, i])*G.val[i + 3N, j + 3N])*
   (-2*G.val[j, j + 2N]*G.val[j + N, i + 3N] - 2*G.val[j, j + N]*G.val[j + 2N, i + 3N] + 
    2*G.val[j, i + 3N]*(G.val[j + N, j + 2N] + G.val[j + 2N, j + N] + G.val[j + 3N, j]) - (-1 + 2*G.val[j, j])*(G.val[j + 3N, i + 3N] - I[j, i])) - 
  (2*G.val[i, i + 2N]*G.val[i + N, j + N] + 2*G.val[i, i + N]*G.val[i + 2N, j + N] - 2*G.val[i, j + N]*(G.val[i + N, i + 2N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) + 
    (-1 + 2*G.val[i, i])*G.val[i + 3N, j + N])*(G.val[j + N, i + 3N]*(-1 + 2*G.val[j + 2N, j + 2N]) + 2*G.val[j, i + 3N]*G.val[j + 2N, j + 3N] - 
    2*G.val[j + 2N, i + 3N]*(G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 3N, j]) + 2*G.val[j + 2N, j]*(G.val[j + 3N, i + 3N] - I[j, i])) - 
  (2*G.val[i, i + 2N]*G.val[i + N, j] + 2*G.val[i, i + N]*G.val[i + 2N, j] - 2*G.val[i, j]*(G.val[i + N, i + 2N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) + 
    (-1 + 2*G.val[i, i])*G.val[i + 3N, j])*(2*G.val[j + 2N, i + 3N]*G.val[j + 3N, j + N] + 2*G.val[j + N, i + 3N]*G.val[j + 3N, j + 2N] + 
    G.val[j, i + 3N]*(-1 + 2*G.val[j + 3N, j + 3N]) - 2*(G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 2N, j + N])*(G.val[j + 3N, i + 3N] - I[j, i])) + 
  4*(G.val[i + N, j + N]*(-1/2 + G.val[i + 2N, i + 2N]) + G.val[i, j + N]*G.val[i + 2N, i + 3N] - 
    G.val[i + 2N, j + N]*(G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 3N, i]) + G.val[i + 2N, i]*G.val[i + 3N, j + N])*
   (G.val[j + N, i + N]*(1/2 - G.val[j + 2N, j + 2N]) - G.val[j, i + N]*G.val[j + 2N, j + 3N] + 
    G.val[j + 2N, i + N]*(G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 3N, j]) - G.val[j + 2N, j]*G.val[j + 3N, i + N] - I[j, i]/2 + 
    G.val[j + 2N, j + 2N]*I[j, i]) + 4*(-((G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 2N, i + N])*G.val[i + 3N, j]) + G.val[i + 2N, j]*G.val[i + 3N, i + N] + 
    G.val[i + N, j]*G.val[i + 3N, i + 2N] + G.val[i, j]*(-1/2 + G.val[i + 3N, i + 3N]))*
   ((G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 2N, j + N])*G.val[j + 3N, i] - G.val[j + 2N, i]*G.val[j + 3N, j + N] - G.val[j + N, i]*G.val[j + 3N, j + 2N] + 
    G.val[j, i]*(1/2 - G.val[j + 3N, j + 3N]) - I[j, i]/2 + G.val[j + 3N, j + 3N]*I[j, i]) + 
  (2*G.val[i + 2N, j + 2N]*G.val[i + 3N, i + N] + 2*G.val[i + N, j + 2N]*G.val[i + 3N, i + 2N] - 2*(G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 2N, i + N])*
     G.val[i + 3N, j + 2N] + G.val[i, j + 2N]*(-1 + 2*G.val[i + 3N, i + 3N]))*((1 - 2*G.val[j + N, j + N])*G.val[j + 2N, i] - 2*G.val[j + N, j]*G.val[j + 3N, i] + 
    2*G.val[j + N, i]*(G.val[j, j + 3N] + G.val[j + 2N, j + N] + G.val[j + 3N, j]) + 2*G.val[j + N, j + 3N]*(-G.val[j, i] + I[j, i])) + 
  4*(G.val[i + 2N, j + N]*G.val[i + 3N, i + N] - (G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 2N, i + N])*G.val[i + 3N, j + N] + 
    G.val[i + N, j + N]*G.val[i + 3N, i + 2N] + G.val[i, j + N]*(-1/2 + G.val[i + 3N, i + 3N]))*(G.val[j + N, i]*(1/2 - G.val[j + 2N, j + 2N]) - 
    G.val[j + 2N, j]*G.val[j + 3N, i] + G.val[j + 2N, i]*(G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 3N, j]) + G.val[j + 2N, j + 3N]*(-G.val[j, i] + I[j, i])) - 
  (2*G.val[i + 2N, j + 3N]*G.val[i + 3N, i + N] + 2*G.val[i + N, j + 3N]*G.val[i + 3N, i + 2N] + G.val[i, j + 3N]*(-1 + 2*G.val[i + 3N, i + 3N]) - 
    2*(G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 2N, i + N])*G.val[i + 3N, j + 3N])*(2*G.val[j, j + 2N]*G.val[j + N, i] + 2*G.val[j, j + N]*G.val[j + 2N, i] + 
    (-1 + 2*G.val[j, j])*G.val[j + 3N, i] + 2*(G.val[j + N, j + 2N] + G.val[j + 2N, j + N] + G.val[j + 3N, j])*(-G.val[j, i] + I[j, i])) + 
  (G.val[i + N, j + 3N]*(-1 + 2*G.val[i + 2N, i + 2N]) + 2*G.val[i, j + 3N]*G.val[i + 2N, i + 3N] - 
    2*G.val[i + 2N, j + 3N]*(G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 3N, i]) + 2*G.val[i + 2N, i]*G.val[i + 3N, j + 3N])*
   (-2*G.val[j, j + N]*G.val[j + 2N, i + N] + 2*G.val[j, i + N]*(G.val[j + N, j + 2N] + G.val[j + 2N, j + N] + G.val[j + 3N, j]) + (1 - 2*G.val[j, j])*G.val[j + 3N, i + N] + 
    2*G.val[j, j + 2N]*(-G.val[j + N, i + N] + I[j, i])) + 4*(G.val[i + N, j]*(-1/2 + G.val[i + 2N, i + 2N]) + G.val[i, j]*G.val[i + 2N, i + 3N] - 
    G.val[i + 2N, j]*(G.val[i, i + 3N] + G.val[i + N, i + 2N] + G.val[i + 3N, i]) + G.val[i + 2N, i]*G.val[i + 3N, j])*
   ((G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 2N, j + N])*G.val[j + 3N, i + N] - G.val[j + 2N, i + N]*G.val[j + 3N, j + N] + 
    G.val[j, i + N]*(1/2 - G.val[j + 3N, j + 3N]) + G.val[j + 3N, j + 2N]*(-G.val[j + N, i + N] + I[j, i])) + 
  4*(G.val[i, j]*G.val[i + N, i + 3N] + (-1/2 + G.val[i + N, i + N])*G.val[i + 2N, j] - G.val[i + N, j]*(G.val[i, i + 3N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) + 
    G.val[i + N, i]*G.val[i + 3N, j])*((G.val[j, j + 3N] + G.val[j + N, j + 2N] + G.val[j + 2N, j + N])*G.val[j + 3N, i + 2N] - 
    G.val[j + N, i + 2N]*G.val[j + 3N, j + 2N] + G.val[j, i + 2N]*(1/2 - G.val[j + 3N, j + 3N]) + G.val[j + 3N, j + N]*(-G.val[j + 2N, i + 2N] + I[j, i])) - 
  (-2*G.val[i, i + 2N]*G.val[i + N, j + 2N] - 2*G.val[i, i + N]*G.val[i + 2N, j + 2N] + 
    2*G.val[i, j + 2N]*(G.val[i + N, i + 2N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) + (1 - 2*G.val[i, i])*G.val[i + 3N, j + 2N])*
   (-2*G.val[j, i + 3N]*G.val[j + N, j + 3N] + (1 - 2*G.val[j + N, j + N])*G.val[j + 2N, i + 3N] + 
    2*G.val[j + N, i + 3N]*(G.val[j, j + 3N] + G.val[j + 2N, j + N] + G.val[j + 3N, j]) + 2*G.val[j + N, j]*(-G.val[j + 3N, i + 3N] + I[j, i])) + 
  4*(G.val[i, j + 2N]*G.val[i + N, i + 3N] + (-1/2 + G.val[i + N, i + N])*G.val[i + 2N, j + 2N] - 
    G.val[i + N, j + 2N]*(G.val[i, i + 3N] + G.val[i + 2N, i + N] + G.val[i + 3N, i]) + G.val[i + N, i]*G.val[i + 3N, j + 2N])*
   (G.val[j + N, i + 2N]*(G.val[j, j + 3N] + G.val[j + 2N, j + N] + G.val[j + 3N, j]) + 
    (-2*G.val[j, i + 2N]*G.val[j + N, j + 3N] + (1 - 2*G.val[j + N, j + N])*G.val[j + 2N, i + 2N] - 2*G.val[j + N, j]*G.val[j + 3N, i + 2N] + 
      (-1 + 2*G.val[j + N, j + N])*I[j, i])/2))/4; 

    return h4a + h4b +4h4c 
end

@inline Base.@propagate_inbounds function heat_cap_h4_kernel(
    mc, model::TwoBandModel,  ij::NTuple{2}, packed_greens::GreensMatrix, 
    flv, ::Discrete_MBF1_X_symm)
    i, j = ij   
	G = packed_greens
    N = length(lattice(model))

  h4a=4*(abs2(G.val[i, i + N]) + abs2(G.val[i + N, i]) + 
  real(G.val[i, i] + 2*conj(G.val[i + N, i])*G.val[i, i + N] + 2*G.val[i, i + N]*G.val[i + N, i] + G.val[i + N, i + N] - 2*G.val[i, i]*G.val[i + N, i + N]))*
 (abs2(G.val[j, j + N]) + abs2(G.val[j + N, j]) + 
  real(G.val[j, j] + 2*conj(G.val[j + N, j])*G.val[j, j + N] + 2*G.val[j, j + N]*G.val[j + N, j] + G.val[j + N, j + N] - 2*G.val[j, j]*G.val[j + N, j + N]));


  h4b=4*(abs2(G.val[i, j + N]*G.val[j, i + N] + G.val[i + N, j]*G.val[j + N, i] + G.val[i + N, j + N]*(G.val[j, i] - I[j, i]) + G.val[i, j]*(G.val[j + N, i + N] - I[j, i])) + 
  2*real((G.val[i, j + N]*G.val[i + N, j] - G.val[i, j]*G.val[i + N, j + N])*(G.val[j, i + N]*G.val[j + N, i] + (G.val[j, i] - I[j, i])*(-G.val[j + N, i + N] + I[j, i]))));
  

  h4c=real(((conj(G.val[i, j])*(-1 + 2*conj(G.val[i + N, i + N])) - 2*conj(G.val[i + N, j])*(G.val[i + N, i] + 2*real(G.val[i, i + N])))*
  (conj(G.val[j, i])*(1 - 2*conj(G.val[j + N, j + N])) + (-1 + 2*conj(G.val[j + N, j + N]))*I[j, i] + 
   2*conj(G.val[j + N, i])*(G.val[j + N, j] + 2*real(G.val[j, j + N]))) + (conj(G.val[j, i + N])*(-1 + 2*conj(G.val[j + N, j + N])) - 
   2*(conj(G.val[j + N, i + N]) - I[j, i])*(G.val[j + N, j] + 2*real(G.val[j, j + N])))*((1 - 2*conj(G.val[i, i]))*conj(G.val[i + N, j]) + 
   2*conj(G.val[i, j])*(G.val[i, i + N] + 2*real(G.val[i + N, i]))) + ((-1 + 2*conj(G.val[i, i]))*conj(G.val[i + N, j + N]) - 
   2*conj(G.val[i, j + N])*(G.val[i, i + N] + 2*real(G.val[i + N, i])))*(-((-1 + 2*conj(G.val[j, j]))*(conj(G.val[j + N, i + N]) - I[j, i])) + 
   2*conj(G.val[j, i + N])*(G.val[j, j + N] + 2*real(G.val[j + N, j]))) - (conj(G.val[i, j + N])*(-1 + 2*conj(G.val[i + N, i + N])) - 
   2*conj(G.val[i + N, j + N])*(G.val[i + N, i] + 2*real(G.val[i, i + N])))*((-1 + 2*conj(G.val[j, j]))*conj(G.val[j + N, i]) - 
   2*(conj(G.val[j, i]) - I[j, i])*(G.val[j, j + N] + 2*real(G.val[j + N, j]))))/2)


  return h4a + h4b +4h4c 
end






