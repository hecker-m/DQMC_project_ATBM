###########################
### lattice iterator for B1 charge density susceptibility
###########################
"""
    EachSitePair_B1()

Creates an iterator template which returns every pair of sites `(s1, s2)` with 
`s1, s2 ∈ 1:Nsites`.
"""
struct EachSitePair_B1 <: DirectLatticeIterator end
EachSitePair_B1(::MonteCarloFlavor) = EachSitePair_B1()
output_size(::EachSitePair_B1, l::Lattice) = (1, )

_length(::EachSitePair_B1, l::Lattice) = 1





function apply!(
    temp::Array, iter::EachSitePair_B1, measurement, mc::DQMC, 
    packed_greens, weight = 1.0
    )
    @timeit_debug "apply!(::EachSitePair_B1, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        srcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L= lat.Ls[1]
                  @inbounds @fastmath for i in eachindex(lat)
                    iPa1=srcdir2trg[i, 2]           # i + a1
                    iMa1=srcdir2trg[i, 1 + L-1]     # i - a1
                    iPa2=srcdir2trg[i, 1 + L]       # i + a2
                    iMa2=srcdir2trg[i,1+L*L-L]  # i - a2

                          @simd for k in eachindex(lat)
                                kPa1=srcdir2trg[k, 2]           # k + a1
                                kMa1=srcdir2trg[k, 1 + L-1]     # k - a1
                                kPa2=srcdir2trg[k, 1 + L]       # k + a2
                                kMa2=srcdir2trg[k, 1 + L*L-L]  # k - a2
                                   
                                temp[1] += 0.25 * weight *measurement.kernel(mc, mc.model, 
                                (i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2), packed_greens, 4)
                          end
                  end
    end
    return 
end


@inline function finalize_temp!(::EachSitePair_B1, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^4
end
@inline function commit!(::EachSitePair_B1, m) 
    push!(m.observable, m.temp[1])
end

###########################
### B1_charge_density susceptibility 
###########################
function B1_charge_density_measurement(
    dqmc::DQMC, model::Model,  greens_iterator; 
        lattice_iterator = EachSitePair_B1(), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = full_B1_charge_density_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end

"""
B1_charge_density_correlation(mc, model; kwargs...)

Generates an equal-time B1_charge_density correlation measurement. Note that the result needs to be added to the simulation 
via `mc[:name] = result`.

## Optional Keyword Arguments
- kwargs from `DQMCMeasurement`
"""
B1_charge_density_correlation(args...; kwargs...) = B1_charge_density_measurement(args..., Greens(); kwargs...)

"""
B1_charge_density_susceptibility(mc, model; kwargs...)

Generates an time-integrated B1_charge_density susceptibility measurement. Note that the result needs to be added to the 
simulation via `mc[:name] = result`.

## Optional Keyword Arguments
- kwargs from `DQMCMeasurement`
"""
B1_charge_density_susceptibility(mc, args...; kwargs...) = B1_charge_density_measurement(mc, args..., TimeIntegral(mc); kwargs...)


###########################
### B1_charge_density susceptibility kernel
###########################


"""
Calculates the B1_charge_density susceptibility kernel 
"""
@inline Base.@propagate_inbounds function full_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10}, G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_B1_charge_density_kernel(mc, model, isks, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_B1_charge_density_kernel(mc::DQMC, model::TwoBandModel, 
    isks::NTuple{10},  G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_B1_charge_density_kernel(mc, model, isks, (G, G, G, G), flv, field)
end



@inline Base.@propagate_inbounds function full_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField  )
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 


    return -(G0l.val[kPa1, i]*Gl0.val[iMa1, k]) + G0l.val[kPa2, i]*Gl0.val[iMa1, k] + G0l.val[kPa1, i]*Gl0.val[iMa2, k] - G0l.val[kPa2, i]*Gl0.val[iMa2, k] - G0l.val[kPa1, i]*Gl0.val[iPa1, k] + 
    G0l.val[kPa2, i]*Gl0.val[iPa1, k] + G0l.val[kMa2, i]*(Gl0.val[iMa1, k] - Gl0.val[iMa2, k] + Gl0.val[iPa1, k] - Gl0.val[iPa2, k]) + G0l.val[kPa1, i]*Gl0.val[iPa2, k] - G0l.val[kPa2, i]*Gl0.val[iPa2, k] + 
    G0l.val[kMa1, i]*(-Gl0.val[iMa1, k] + Gl0.val[iMa2, k] - Gl0.val[iPa1, k] + Gl0.val[iPa2, k]) - (G0l.val[kMa1 + N, i] - G0l.val[kMa2 + N, i] + G0l.val[kPa1 + N, i] - G0l.val[kPa2 + N, i])*
     (Gl0.val[iMa1, k + N] - Gl0.val[iMa2, k + N] + Gl0.val[iPa1, k + N] - Gl0.val[iPa2, k + N]) - (G0l.val[kMa1 + 2N, i] - G0l.val[kMa2 + 2N, i] + G0l.val[kPa1 + 2N, i] - G0l.val[kPa2 + 2N, i])*
     (Gl0.val[iMa1, k + 2N] - Gl0.val[iMa2, k + 2N] + Gl0.val[iPa1, k + 2N] - Gl0.val[iPa2, k + 2N]) - 
    (G0l.val[kMa1 + 3N, i] - G0l.val[kMa2 + 3N, i] + G0l.val[kPa1 + 3N, i] - G0l.val[kPa2 + 3N, i])*(Gl0.val[iMa1, k + 3N] - Gl0.val[iMa2, k + 3N] + Gl0.val[iPa1, k + 3N] - 
      Gl0.val[iPa2, k + 3N]) - (G0l.val[kMa1, i + N] - G0l.val[kMa2, i + N] + G0l.val[kPa1, i + N] - G0l.val[kPa2, i + N])*(Gl0.val[iMa1 + N, k] - Gl0.val[iMa2 + N, k] + Gl0.val[iPa1 + N, k] - 
      Gl0.val[iPa2 + N, k]) - (G0l.val[kMa1 + N, i + N] - G0l.val[kMa2 + N, i + N] + G0l.val[kPa1 + N, i + N] - G0l.val[kPa2 + N, i + N])*
     (Gl0.val[iMa1 + N, k + N] - Gl0.val[iMa2 + N, k + N] + Gl0.val[iPa1 + N, k + N] - Gl0.val[iPa2 + N, k + N]) - 
    (G0l.val[kMa1 + 2N, i + N] - G0l.val[kMa2 + 2N, i + N] + G0l.val[kPa1 + 2N, i + N] - G0l.val[kPa2 + 2N, i + N])*
     (Gl0.val[iMa1 + N, k + 2N] - Gl0.val[iMa2 + N, k + 2N] + Gl0.val[iPa1 + N, k + 2N] - Gl0.val[iPa2 + N, k + 2N]) - 
    (G0l.val[kMa1 + 3N, i + N] - G0l.val[kMa2 + 3N, i + N] + G0l.val[kPa1 + 3N, i + N] - G0l.val[kPa2 + 3N, i + N])*
     (Gl0.val[iMa1 + N, k + 3N] - Gl0.val[iMa2 + N, k + 3N] + Gl0.val[iPa1 + N, k + 3N] - Gl0.val[iPa2 + N, k + 3N]) - 
    (G0l.val[kMa1, i + 2N] - G0l.val[kMa2, i + 2N] + G0l.val[kPa1, i + 2N] - G0l.val[kPa2, i + 2N])*(Gl0.val[iMa1 + 2N, k] - Gl0.val[iMa2 + 2N, k] + Gl0.val[iPa1 + 2N, k] - 
      Gl0.val[iPa2 + 2N, k]) - (G0l.val[kMa1 + N, i + 2N] - G0l.val[kMa2 + N, i + 2N] + G0l.val[kPa1 + N, i + 2N] - G0l.val[kPa2 + N, i + 2N])*
     (Gl0.val[iMa1 + 2N, k + N] - Gl0.val[iMa2 + 2N, k + N] + Gl0.val[iPa1 + 2N, k + N] - Gl0.val[iPa2 + 2N, k + N]) - 
    (G0l.val[kMa1 + 2N, i + 2N] - G0l.val[kMa2 + 2N, i + 2N] + G0l.val[kPa1 + 2N, i + 2N] - G0l.val[kPa2 + 2N, i + 2N])*
     (Gl0.val[iMa1 + 2N, k + 2N] - Gl0.val[iMa2 + 2N, k + 2N] + Gl0.val[iPa1 + 2N, k + 2N] - Gl0.val[iPa2 + 2N, k + 2N]) - 
    (G0l.val[kMa1 + 3N, i + 2N] - G0l.val[kMa2 + 3N, i + 2N] + G0l.val[kPa1 + 3N, i + 2N] - G0l.val[kPa2 + 3N, i + 2N])*
     (Gl0.val[iMa1 + 2N, k + 3N] - Gl0.val[iMa2 + 2N, k + 3N] + Gl0.val[iPa1 + 2N, k + 3N] - Gl0.val[iPa2 + 2N, k + 3N]) - 
    (G0l.val[kMa1, i + 3N] - G0l.val[kMa2, i + 3N] + G0l.val[kPa1, i + 3N] - G0l.val[kPa2, i + 3N])*(Gl0.val[iMa1 + 3N, k] - Gl0.val[iMa2 + 3N, k] + Gl0.val[iPa1 + 3N, k] - 
      Gl0.val[iPa2 + 3N, k]) - (G0l.val[kMa1 + N, i + 3N] - G0l.val[kMa2 + N, i + 3N] + G0l.val[kPa1 + N, i + 3N] - G0l.val[kPa2 + N, i + 3N])*
     (Gl0.val[iMa1 + 3N, k + N] - Gl0.val[iMa2 + 3N, k + N] + Gl0.val[iPa1 + 3N, k + N] - Gl0.val[iPa2 + 3N, k + N]) - 
    (G0l.val[kMa1 + 2N, i + 3N] - G0l.val[kMa2 + 2N, i + 3N] + G0l.val[kPa1 + 2N, i + 3N] - G0l.val[kPa2 + 2N, i + 3N])*
     (Gl0.val[iMa1 + 3N, k + 2N] - Gl0.val[iMa2 + 3N, k + 2N] + Gl0.val[iPa1 + 3N, k + 2N] - Gl0.val[iPa2 + 3N, k + 2N]) - 
    (G0l.val[kMa1 + 3N, i + 3N] - G0l.val[kMa2 + 3N, i + 3N] + G0l.val[kPa1 + 3N, i + 3N] - G0l.val[kPa2 + 3N, i + 3N])*
     (Gl0.val[iMa1 + 3N, k + 3N] - Gl0.val[iMa2 + 3N, k + 3N] + Gl0.val[iPa1 + 3N, k + 3N] - Gl0.val[iPa2 + 3N, k + 3N]) + 
    (G00.val[kMa1, k] - G00.val[kMa2, k] + G00.val[kPa1, k] - G00.val[kPa2, k] + G00.val[kMa1 + N, k + N] - G00.val[kMa2 + N, k + N] + G00.val[kPa1 + N, k + N] - G00.val[kPa2 + N, k + N] + 
      G00.val[kMa1 + 2N, k + 2N] - G00.val[kMa2 + 2N, k + 2N] + G00.val[kPa1 + 2N, k + 2N] - G00.val[kPa2 + 2N, k + 2N] + G00.val[kMa1 + 3N, k + 3N] - G00.val[kMa2 + 3N, k + 3N] + 
      G00.val[kPa1 + 3N, k + 3N] - G00.val[kPa2 + 3N, k + 3N])*(Gll.val[iMa1, i] - Gll.val[iMa2, i] + Gll.val[iPa1, i] - Gll.val[iPa2, i] + Gll.val[iMa1 + N, i + N] - Gll.val[iMa2 + N, i + N] + 
      Gll.val[iPa1 + N, i + N] - Gll.val[iPa2 + N, i + N] + Gll.val[iMa1 + 2N, i + 2N] - Gll.val[iMa2 + 2N, i + 2N] + Gll.val[iPa1 + 2N, i + 2N] - Gll.val[iPa2 + 2N, i + 2N] + 
      Gll.val[iMa1 + 3N, i + 3N] - Gll.val[iMa2 + 3N, i + 3N] + Gll.val[iPa1 + 3N, i + 3N] - Gll.val[iPa2 + 3N, i + 3N]) + id*Gl0.val[iMa1, k]*I[kMa1, i] - 
    id*Gl0.val[iMa2, k]*I[kMa1, i] + id*Gl0.val[iPa1, k]*I[kMa1, i] - id*Gl0.val[iPa2, k]*I[kMa1, i] + id*Gl0.val[iMa1 + N, k + N]*I[kMa1, i] - id*Gl0.val[iMa2 + N, k + N]*I[kMa1, i] + 
    id*Gl0.val[iPa1 + N, k + N]*I[kMa1, i] - id*Gl0.val[iPa2 + N, k + N]*I[kMa1, i] + id*Gl0.val[iMa1 + 2N, k + 2N]*I[kMa1, i] - id*Gl0.val[iMa2 + 2N, k + 2N]*I[kMa1, i] + 
    id*Gl0.val[iPa1 + 2N, k + 2N]*I[kMa1, i] - id*Gl0.val[iPa2 + 2N, k + 2N]*I[kMa1, i] + id*Gl0.val[iMa1 + 3N, k + 3N]*I[kMa1, i] - 
    id*Gl0.val[iMa2 + 3N, k + 3N]*I[kMa1, i] + id*Gl0.val[iPa1 + 3N, k + 3N]*I[kMa1, i] - id*Gl0.val[iPa2 + 3N, k + 3N]*I[kMa1, i] - id*Gl0.val[iMa1, k]*I[kMa2, i] + 
    id*Gl0.val[iMa2, k]*I[kMa2, i] - id*Gl0.val[iPa1, k]*I[kMa2, i] + id*Gl0.val[iPa2, k]*I[kMa2, i] - id*Gl0.val[iMa1 + N, k + N]*I[kMa2, i] + id*Gl0.val[iMa2 + N, k + N]*I[kMa2, i] - 
    id*Gl0.val[iPa1 + N, k + N]*I[kMa2, i] + id*Gl0.val[iPa2 + N, k + N]*I[kMa2, i] - id*Gl0.val[iMa1 + 2N, k + 2N]*I[kMa2, i] + id*Gl0.val[iMa2 + 2N, k + 2N]*I[kMa2, i] - 
    id*Gl0.val[iPa1 + 2N, k + 2N]*I[kMa2, i] + id*Gl0.val[iPa2 + 2N, k + 2N]*I[kMa2, i] - id*Gl0.val[iMa1 + 3N, k + 3N]*I[kMa2, i] + 
    id*Gl0.val[iMa2 + 3N, k + 3N]*I[kMa2, i] - id*Gl0.val[iPa1 + 3N, k + 3N]*I[kMa2, i] + id*Gl0.val[iPa2 + 3N, k + 3N]*I[kMa2, i] + id*Gl0.val[iMa1, k]*I[kPa1, i] - 
    id*Gl0.val[iMa2, k]*I[kPa1, i] + id*Gl0.val[iPa1, k]*I[kPa1, i] - id*Gl0.val[iPa2, k]*I[kPa1, i] + id*Gl0.val[iMa1 + N, k + N]*I[kPa1, i] - id*Gl0.val[iMa2 + N, k + N]*I[kPa1, i] + 
    id*Gl0.val[iPa1 + N, k + N]*I[kPa1, i] - id*Gl0.val[iPa2 + N, k + N]*I[kPa1, i] + id*Gl0.val[iMa1 + 2N, k + 2N]*I[kPa1, i] - id*Gl0.val[iMa2 + 2N, k + 2N]*I[kPa1, i] + 
    id*Gl0.val[iPa1 + 2N, k + 2N]*I[kPa1, i] - id*Gl0.val[iPa2 + 2N, k + 2N]*I[kPa1, i] + id*Gl0.val[iMa1 + 3N, k + 3N]*I[kPa1, i] - 
    id*Gl0.val[iMa2 + 3N, k + 3N]*I[kPa1, i] + id*Gl0.val[iPa1 + 3N, k + 3N]*I[kPa1, i] - id*Gl0.val[iPa2 + 3N, k + 3N]*I[kPa1, i] + 
    id*(-Gl0.val[iMa1, k] + Gl0.val[iMa2, k] - Gl0.val[iPa1, k] + Gl0.val[iPa2, k] - Gl0.val[iMa1 + N, k + N] + Gl0.val[iMa2 + N, k + N] - Gl0.val[iPa1 + N, k + N] + Gl0.val[iPa2 + N, k + N] - 
      Gl0.val[iMa1 + 2N, k + 2N] + Gl0.val[iMa2 + 2N, k + 2N] - Gl0.val[iPa1 + 2N, k + 2N] + Gl0.val[iPa2 + 2N, k + 2N] - Gl0.val[iMa1 + 3N, k + 3N] + Gl0.val[iMa2 + 3N, k + 3N] - 
      Gl0.val[iPa1 + 3N, k + 3N] + Gl0.val[iPa2 + 3N, k + 3N])*I[kPa2, i]
end

@inline Base.@propagate_inbounds function full_B1_charge_density_kernel(
        mc, model::TwoBandModel, isks::NTuple{10}, packed_greens::_GM4{<: Matrix}, 
        flv, ::Union{Discrete_MBF1_X_symm, Discrete_MBF2_symm})
    i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2 = isks   # i, i+a1, i-a1, etc.
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 


    return real(2*G0l.val[kMa2, i]*(Gl0.val[iMa1, k] - Gl0.val[iMa2, k] + Gl0.val[iPa1, k] - Gl0.val[iPa2, k]) + 
    2*G0l.val[kMa1, i]*(-Gl0.val[iMa1, k] + Gl0.val[iMa2, k] - Gl0.val[iPa1, k] + Gl0.val[iPa2, k]) + 
    (G0l.val[kMa1 + N, i] - G0l.val[kMa2 + N, i] + G0l.val[kPa1 + N, i] - G0l.val[kPa2 + N, i])*(Gl0.val[iMa1, k + N] - Gl0.val[iMa2, k + N] + Gl0.val[iPa1, k + N] - 
    Gl0.val[iPa2, k + N]) + (G0l.val[kMa1, i + N] - G0l.val[kMa2, i + N] + G0l.val[kPa1, i + N] - G0l.val[kPa2, i + N])*
    (Gl0.val[iMa1 + N, k] - Gl0.val[iMa2 + N, k] + Gl0.val[iPa1 + N, k] - Gl0.val[iPa2 + N, k]) + 
    (G0l.val[kMa1 + N, i + N] - G0l.val[kMa2 + N, i + N] + G0l.val[kPa1 + N, i + N] - G0l.val[kPa2 + N, i + N])*
    (Gl0.val[iMa1 + N, k + N] - Gl0.val[iMa2 + N, k + N] + Gl0.val[iPa1 + N, k + N] - Gl0.val[iPa2 + N, k + N]) + 
    2*(G0l.val[kPa2, i]*(Gl0.val[iMa1, k] - Gl0.val[iMa2, k] + Gl0.val[iPa1, k] - Gl0.val[iPa2, k]) + 
    G0l.val[kPa1, i]*(-Gl0.val[iMa1, k] + Gl0.val[iMa2, k] - Gl0.val[iPa1, k] + Gl0.val[iPa2, k]) + 
    id*(Gl0.val[iMa1, k] - Gl0.val[iMa2, k] + Gl0.val[iPa1, k] - Gl0.val[iPa2, k] + Gl0.val[iMa1 + N, k + N] - Gl0.val[iMa2 + N, k + N] + Gl0.val[iPa1 + N, k + N] - 
        Gl0.val[iPa2 + N, k + N])*(I[kMa1, i] - I[kMa2, i] + I[kPa1, i] - I[kPa2, i])))+
        4*real(G00.val[kMa1, k] - G00.val[kMa2, k] + G00.val[kPa1, k] - G00.val[kPa2, k] + G00.val[kMa1 + N, k + N] - G00.val[kMa2 + N, k + N] + G00.val[kPa1 + N, k + N] - 
        G00.val[kPa2 + N, k + N])*real(Gll.val[iMa1, i] - Gll.val[iMa2, i] + Gll.val[iPa1, i] - Gll.val[iPa2, i] + Gll.val[iMa1 + N, i + N] - Gll.val[iMa2 + N, i + N] + 
        Gll.val[iPa1 + N, i + N] - Gll.val[iPa2 + N, i + N]) 
end


