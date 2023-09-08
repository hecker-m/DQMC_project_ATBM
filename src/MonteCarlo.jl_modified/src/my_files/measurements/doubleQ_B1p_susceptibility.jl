###########################
### lattice iterator for B1`_Q1Q2 susceptibility
###########################

"""
    EachSitePairByDistance()

Creates an iterator template which returns triplets 
`(direction index, source, target)` sorted by distance. The `direction index` 
identifies each unique direction `position(target) - position(source)`.

Requires `lattice` to implement `positions` and `lattice_vectors`.
"""
struct EachDoubleSitePairByDistance_B1p_Q1Q2 <: DeferredLatticeIterator end
EachDoubleSitePairByDistance_B1p_Q1Q2(::MonteCarloFlavor) = EachDoubleSitePairByDistance_B1p_Q1Q2()

function output_size(::EachDoubleSitePairByDistance_B1p_Q1Q2, l::Lattice)
    return (1,)
end



function apply!(
    temp::Array, iter::EachDoubleSitePairByDistance_B1p_Q1Q2, measurement, mc::DQMC, 
    packed_greens, weight = 1.0
    )
    @timeit_debug "apply!(::EachDoubleSitePairByDistance_B1p_Q1Q2, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L= lat.Ls[1]

        k=div(N,2)
         @inbounds @fastmath for ix in 1:L, iy in 1:L
              if isodd(ix+iy)
                  Δi =(-1)^ix -(-1)^iy  
                  #note that i and j are treated as directions [so is p], 
                  #i.e. in principle, we should subtract (ix,iy) = (ix, iy)- (1, 1). Easy to see, that these subtractions cancel.
                  i=ix+(iy-1)*L                   
                  kPi = Bsrcdir2trg[k, i]
                  for jx in 1:L, jy in 1:L
                      if isodd(jx+jy)
                          Δj=(-1)^jx -(-1)^jy
                          j=jx+(jy-1)*L
                          tcc=zero(ComplexF64);
                          @simd for p in eachindex(lat)
                                    kPp=Bsrcdir2trg[k, p]        #l=k+p
                                    kPpPj = Bsrcdir2trg[kPp, j]  #l+j=k+p+j
                                    py, px=fldmod1(p, L)                                    
                                  tcc +=(-1)^(px+py) *measurement.kernel(mc, mc.model, (k, kPp, kPi, kPpPj), packed_greens, 4)
                          end
                          temp[1] += weight * Δi*Δj*tcc*N
                      end
                  end
              end
          end
    end
    return 
end


@inline function finalize_temp!(::EachDoubleSitePairByDistance_B1p_Q1Q2, m, mc)
    m.temp ./= (length(lattice(mc))^4)  # *1/N^4
end
@inline function commit!(::EachDoubleSitePairByDistance_B1p_Q1Q2, m) 
    push!(m.observable, m.temp[1])
end

###########################
### B1p_Q1Q2 double-Q susceptibility 
###########################
function B1p_Q1Q2_measurement(
    dqmc::DQMC, model::Model,  greens_iterator; 
        lattice_iterator = EachDoubleSitePairByDistance_B1p_Q1Q2(), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = full_B1p_Q1Q2_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end

"""
B1`_Q1Q2_correlation(mc, model; kwargs...)

Generates an equal-time B1`_Q1Q2 double-Q correlation measurement. Note that the result needs to be added to the simulation 
via `mc[:name] = result`.

## Optional Keyword Arguments
- kwargs from `DQMCMeasurement`
"""
B1p_Q1Q2_correlation(args...; kwargs...) = B1p_Q1Q2_measurement(args..., Greens(); kwargs...)

"""
B1`_Q1Q2_susceptibility(mc, model; kwargs...)

Generates an time-integrated B1`_Q1Q2 double-Q susceptibility measurement. Note that the result needs to be added to the 
simulation via `mc[:name] = result`.

## Optional Keyword Arguments
- kwargs from `DQMCMeasurement`
"""
B1p_Q1Q2_susceptibility(mc, args...; kwargs...) = B1p_Q1Q2_measurement(mc, args..., TimeIntegral(mc); kwargs...)


###########################
### B1` double-Q susceptibility kernel
###########################


"""
Calculates the B1` double-Q` susceptibility kernel 
"""
@inline Base.@propagate_inbounds function full_B1p_Q1Q2_kernel(mc::DQMC, model::TwoBandModel, klkPlP::NTuple{4}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_B1p_Q1Q2_kernel(mc, model, klkPlP, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_B1p_Q1Q2_kernel(mc::DQMC, model::TwoBandModel, klkPlP::NTuple{4}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_B1p_Q1Q2_kernel(mc, model, klkPlP, (G, G, G, G), flv, field)
end



@inline Base.@propagate_inbounds function full_B1p_Q1Q2_kernel(
        mc, model::TwoBandModel, klkPlP::NTuple{4}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
    )
    k, l, kPi, lPj = klkPlP   #k, l, k+i, l+j
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[G0l.k, G0l.l] 

    tcoscos1=-((-((Gll.val[k, kPi + 2N] + Gll.val[k + N, kPi + 3N])*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N])) + (Gll.val[k, kPi + 3N] + Gll.val[k + N, kPi + 2N])*
    (Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + (Gll.val[kPi, k] + Gll.val[kPi + N, k + N])*(Gll.val[k + 2N, kPi + 3N] + Gll.val[k + 3N, kPi + 2N]) - 
   (Gll.val[kPi, k + N] + Gll.val[kPi + N, k])*(Gll.val[k + 2N, kPi + 2N] + Gll.val[k + 3N, kPi + 3N]) - (Gll.val[k + 2N, kPi] + Gll.val[k + 3N, kPi + N])*
    (Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) - (Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N])*
    (Gll.val[kPi, kPi + 3N] - Gll.val[kPi + N, kPi + 2N] + Gll.val[kPi + 2N, kPi + N] - Gll.val[kPi + 3N, kPi]) + 
   (Gll.val[k + 2N, kPi + N] + Gll.val[k + 3N, kPi])*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) + 
   (Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k])*(Gll.val[kPi, kPi + 2N] - Gll.val[kPi + N, kPi + 3N] + Gll.val[kPi + 2N, kPi] - 
     Gll.val[kPi + 3N, kPi + N]) - (Gll.val[k, kPi] + Gll.val[k + N, kPi + N])*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + 
   (Gll.val[k, kPi + N] + Gll.val[k + N, kPi])*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]) + 
   2*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k] + Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N])*I[k, kPi])*
  (-((G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*(G00.val[lPj, l + 3N] + G00.val[lPj + N, l + 2N])) + (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*
    (G00.val[lPj, l + 2N] + G00.val[lPj + N, l + 3N]) + (G00.val[lPj, l] + G00.val[lPj + N, l + N])*(G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N]) - 
   (G00.val[lPj, l + N] + G00.val[lPj + N, l])*(G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N]) - (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*
    (G00.val[lPj + 2N, l + N] + G00.val[lPj + 3N, l]) - (G00.val[l, l + 2N] - G00.val[l + N, l + 3N] + G00.val[l + 2N, l] - G00.val[l + 3N, l + N])*
    (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj]) + 
   (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*(G00.val[lPj + 2N, l] + G00.val[lPj + 3N, l + N]) + 
   (G00.val[l, l + 3N] - G00.val[l + N, l + 2N] + G00.val[l + 2N, l + N] - G00.val[l + 3N, l])*(G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - 
     G00.val[lPj + 3N, lPj + N]) - (G00.val[l, lPj] + G00.val[l + N, lPj + N])*(G00.val[lPj + 2N, l + 3N] + G00.val[lPj + 3N, l + 2N]) + 
   (G00.val[l, lPj + N] + G00.val[l + N, lPj])*(G00.val[lPj + 2N, l + 2N] + G00.val[lPj + 3N, l + 3N]) + 
   2*(G00.val[lPj, l + N] + G00.val[lPj + N, l] + G00.val[lPj + 2N, l + 3N] + G00.val[lPj + 3N, l + 2N])*I[l, lPj]));

    tcoscos2=4*((G0l.val[l, kPi + 2N]*G0l.val[lPj, k + 2N] - G0l.val[l, k + 2N]*G0l.val[lPj, kPi + 2N] - G0l.val[l, kPi + 3N]*G0l.val[lPj, k + 3N] + G0l.val[l, k + 3N]*G0l.val[lPj, kPi + 3N] - 
    G0l.val[l + N, kPi + 2N]*G0l.val[lPj + N, k + 2N] + G0l.val[l + N, k + 2N]*G0l.val[lPj + N, kPi + 2N] + G0l.val[l + N, kPi + 3N]*G0l.val[lPj + N, k + 3N] - 
    G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi + 3N])*(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + 2N] - Gl0.val[kPi, lPj + 2N]*Gl0.val[k + N, l + 3N]) - 
  (G0l.val[lPj, kPi + 2N]*G0l.val[l + N, k + 2N] - G0l.val[lPj, k + 2N]*G0l.val[l + N, kPi + 2N] - G0l.val[lPj, kPi + 3N]*G0l.val[l + N, k + 3N] + 
    G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi + 3N])*(Gl0.val[kPi, lPj + 2N]*Gl0.val[k + N, l + 2N] - Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + 2N] - 
    Gl0.val[kPi, lPj + 3N]*Gl0.val[k + N, l + 3N] + Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + 3N]) - 
  (G0l.val[l, k + 3N]*G0l.val[lPj, kPi + 2N] - G0l.val[l, kPi + 2N]*G0l.val[lPj, k + 3N] - G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi + 2N] + 
    G0l.val[l + N, kPi + 2N]*G0l.val[lPj + N, k + 3N])*(Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi, l + 3N] - 
    Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + 2N] + Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + N, l + 3N]) + 
  (G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi + 2N] - G0l.val[lPj, kPi + 2N]*G0l.val[l + N, k + 3N])*(Gl0.val[k, lPj + 2N]*Gl0.val[kPi, l + 2N] - Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + 2N] - 
    Gl0.val[k, lPj + 3N]*Gl0.val[kPi, l + 3N] + Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + 3N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + N, l + 2N] + 
    Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + 2N] + Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + N, l + 3N] - Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + 3N]) + 
  (G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k] - G0l.val[l + 2N, k]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 2N, k + N] + 
    G0l.val[l + 2N, k + N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k] + G0l.val[l + 3N, k]*G0l.val[lPj + 3N, kPi] + 
    G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + N] - G0l.val[l + 3N, k + N]*G0l.val[lPj + 3N, kPi + N])*(Gl0.val[kPi + 2N, l + N]*Gl0.val[k + 3N, lPj] - 
    Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + N]) - (G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k] - G0l.val[lPj + 2N, k]*G0l.val[l + 3N, kPi] - 
    G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + N] + G0l.val[lPj + 2N, k + N]*G0l.val[l + 3N, kPi + N])*(Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l] - 
    Gl0.val[kPi + 2N, l]*Gl0.val[k + 3N, lPj] - Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + N] + Gl0.val[kPi + 2N, l + N]*Gl0.val[k + 3N, lPj + N]) - 
  (G0l.val[l + 2N, k + N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + N] - G0l.val[l + 3N, k + N]*G0l.val[lPj + 3N, kPi] + 
    G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k + N])*(Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + N] - 
    Gl0.val[k + 3N, l + N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + N]) + 
  (G0l.val[lPj + 2N, k + N]*G0l.val[l + 3N, kPi] - G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k + N])*(Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l] - Gl0.val[k + 2N, l]*Gl0.val[kPi + 2N, lPj] - 
    Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + N] + Gl0.val[k + 2N, l + N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l] + 
    Gl0.val[k + 3N, l]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + N] - Gl0.val[k + 3N, l + N]*Gl0.val[kPi + 3N, lPj + N]) + 
  (Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + 2N] + 
    Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] + 
    Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + 2N] - Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(-(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + N]) - 
    G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi] + G0l.val[lPj + 2N, kPi]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[lPj + 3N, k + N]*(G0l.val[l, kPi] - id*I[l, kPi])) - 
  (-(Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj]) + Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + N] + Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 2N] - 
    Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 3N])*(-(G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k]) + G0l.val[l + N, k]*G0l.val[lPj + 2N, kPi] - 
    G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + N] + G0l.val[l, k + N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[lPj + 3N, kPi]*(-G0l.val[l, k] + id*I[l, k]) + 
    G0l.val[lPj + 2N, kPi + N]*(-G0l.val[l + N, k + N] + id*I[l, k]) + G0l.val[lPj + 3N, k]*(G0l.val[l, kPi] - id*I[l, kPi]) + 
    G0l.val[lPj + 2N, k + N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi])) - (Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] + 
    Gl0.val[k + 2N, lPj + N]*Gl0.val[kPi + 2N, l + 2N] - Gl0.val[k + 2N, lPj]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj] + 
    Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + 3N, lPj + N]*Gl0.val[kPi + 3N, l + 2N] + Gl0.val[k + 3N, lPj]*Gl0.val[kPi + 3N, l + 3N])*
   (G0l.val[l, k + N]*G0l.val[lPj + 2N, kPi] + G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k + N] + G0l.val[lPj + 3N, kPi]*(-G0l.val[l + N, k + N] + id*I[l, k]) + 
    G0l.val[lPj + 2N, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi])) + (-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj]) + Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + N] - 
    Gl0.val[kPi + 2N, lPj + N]*Gl0.val[k + 3N, l + 2N] + Gl0.val[kPi + 2N, lPj]*Gl0.val[k + 3N, l + 3N])*(G0l.val[l, kPi + N]*G0l.val[lPj + 2N, k + N] - 
    G0l.val[l, k + N]*G0l.val[lPj + 2N, kPi + N] + G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k] - G0l.val[l + N, k]*G0l.val[lPj + 3N, kPi] + 
    G0l.val[lPj + 2N, kPi]*(G0l.val[l, k] - id*I[l, k]) + G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + N, k + N] - id*I[l, k]) + G0l.val[lPj + 2N, k]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
    G0l.val[lPj + 3N, k + N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi])) + (-(Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l]) + Gl0.val[k, l]*Gl0.val[kPi + 2N, lPj] + 
    Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi + 2N, lPj + N] + Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l] - Gl0.val[k + N, l]*Gl0.val[kPi + 3N, lPj] - 
    Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + N] + Gl0.val[k + N, l + N]*Gl0.val[kPi + 3N, lPj + N])*(-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi]) - 
    G0l.val[lPj + 2N, kPi + N]*G0l.val[l + 3N, k + 2N] + G0l.val[lPj + 2N, kPi]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + 
    G0l.val[l + 3N, kPi + N]*(G0l.val[lPj + 2N, k + 2N] - id*I[lPj, k])) - (Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + N] + 
    Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 2N] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj] + 
    Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + 2N] + Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + 3N])*
   (G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l, k + 2N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi] + 
    G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
    G0l.val[lPj + 3N, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + G0l.val[l, kPi + N]*(G0l.val[lPj + 2N, k + 2N] - id*I[lPj, k]) + 
    G0l.val[l + N, kPi]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k])) - (Gl0.val[k, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + N] - 
    Gl0.val[k + N, l + N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + N])*(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi] - 
    G0l.val[l + 2N, kPi]*G0l.val[lPj + 2N, k + 3N] + G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + N] - G0l.val[l + 3N, kPi + N]*G0l.val[lPj + 3N, k + 2N] + 
    G0l.val[lPj + 2N, kPi + N]*(-G0l.val[l + 2N, k + 2N] + id*I[l, k]) + G0l.val[lPj + 3N, kPi]*(-G0l.val[l + 3N, k + 3N] + id*I[l, k]) + 
    G0l.val[l + 2N, kPi + N]*(G0l.val[lPj + 2N, k + 2N] - id*I[lPj, k]) + G0l.val[l + 3N, kPi]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k])) - 
  (-(Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l]) + Gl0.val[k + N, l]*Gl0.val[kPi + 2N, lPj] + Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + N] - 
    Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj + N] + Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l] - Gl0.val[k, l]*Gl0.val[kPi + 3N, lPj] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + N] + 
    Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj + N])*(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi + N] + G0l.val[lPj + 2N, kPi]*G0l.val[l + 3N, k + 2N] + 
    G0l.val[lPj + 2N, kPi + N]*(-G0l.val[l + 3N, k + 3N] + id*I[l, k]) + G0l.val[l + 3N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k])) - 
  (Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 2N] + 
    Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + N] + 
    Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 2N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi] - 
    G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + N] - G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi] + G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + N] + 
    G0l.val[lPj + 3N, k + 2N]*(G0l.val[l, kPi] - id*I[l, kPi]) + G0l.val[lPj + 2N, k + 3N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) + 
    G0l.val[l + N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  (Gl0.val[k + N, l + N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + N])*
   (-(G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi + N]) + G0l.val[l + 2N, kPi + N]*G0l.val[lPj + 2N, k + 3N] - G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi] + 
    G0l.val[l + 3N, kPi]*G0l.val[lPj + 3N, k + 2N] + G0l.val[lPj + 2N, kPi]*(G0l.val[l + 2N, k + 2N] - id*I[l, k]) + 
    G0l.val[lPj + 3N, kPi + N]*(G0l.val[l + 3N, k + 3N] - id*I[l, k]) + G0l.val[l + 2N, kPi]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
    G0l.val[l + 3N, kPi + N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + (Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + N] + 
    Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 2N, l + 2N] - Gl0.val[k + N, lPj]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj] + 
    Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + N] - Gl0.val[k, lPj + N]*Gl0.val[kPi + 3N, l + 2N] + Gl0.val[k, lPj]*Gl0.val[kPi + 3N, l + 3N])*
   (G0l.val[l, k + 2N]*G0l.val[lPj + 2N, kPi] - G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi + N] + G0l.val[l, kPi + N]*G0l.val[lPj + 2N, k + 3N] - 
    G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi] + G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi + N] + G0l.val[l + N, kPi]*G0l.val[lPj + 3N, k + 2N] + 
    (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) + 
  (Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj] - Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + N] - Gl0.val[k, lPj]*Gl0.val[kPi + 2N, l + 2N] + Gl0.val[k, lPj + N]*Gl0.val[kPi + 2N, l + 3N] - 
    Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj] + Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + N] + Gl0.val[k + N, lPj]*Gl0.val[kPi + 3N, l + 2N] - 
    Gl0.val[k + N, lPj + N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi] - G0l.val[l + N, k + 2N]*G0l.val[lPj + 2N, kPi + N] - 
    G0l.val[l + N, kPi]*G0l.val[lPj + 2N, k + 3N] - G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi] + G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + N] - 
    G0l.val[l, kPi + N]*G0l.val[lPj + 3N, k + 2N] + (-G0l.val[l + N, kPi + N] + id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
    (-G0l.val[l, kPi] + id*I[l, kPi])*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k])) - 
  (Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] + 
    Gl0.val[k + 3N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*(G0l.val[lPj + N, kPi]*(-G0l.val[l + N, k + N] + id*I[l, k]) + G0l.val[lPj, k + N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
    G0l.val[l + N, kPi]*(G0l.val[lPj + N, k + N] - id*I[lPj, k]) + G0l.val[l, k + N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
  (Gl0.val[kPi + 2N, lPj + 2N]*Gl0.val[k + 3N, l + 2N] - Gl0.val[kPi + 2N, l + 2N]*Gl0.val[k + 3N, lPj + 2N] - Gl0.val[kPi + 2N, lPj + 3N]*Gl0.val[k + 3N, l + 3N] + 
    Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + 3N])*(G0l.val[lPj, kPi + N]*(-G0l.val[l + N, k + N] + id*I[l, k]) + 
    G0l.val[lPj, k + N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) + G0l.val[l + N, kPi]*(-G0l.val[lPj, k] + id*I[lPj, k]) + G0l.val[l + N, k]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
  (-(Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N] + Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N] - 
    Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N] + Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N] - Gl0.val[k, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N] - 
    Gl0.val[k, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N] + Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*(-(G0l.val[lPj, k + 2N]*G0l.val[l + N, kPi]) - 
    G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 3N] + G0l.val[lPj, k + 3N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) + G0l.val[l + N, k + 2N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
  (-(Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + Gl0.val[k, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N] + Gl0.val[k, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N] - 
    Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N] + Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N] - Gl0.val[k + N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N] - 
    Gl0.val[k + N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N] + Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*
   (-(G0l.val[lPj, k + 3N]*G0l.val[l + N, kPi]) - G0l.val[lPj, kPi + N]*G0l.val[l + N, k + 2N] + G0l.val[lPj, k + 2N]*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) + 
    G0l.val[l + N, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + (-(Gl0.val[k + 2N, lPj + 2N]*Gl0.val[kPi + 2N, l + 2N]) + 
    Gl0.val[k + 2N, l + 2N]*Gl0.val[kPi + 2N, lPj + 2N] + Gl0.val[k + 2N, lPj + 3N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + 2N, l + 3N]*Gl0.val[kPi + 2N, lPj + 3N] + 
    Gl0.val[k + 3N, lPj + 2N]*Gl0.val[kPi + 3N, l + 2N] - Gl0.val[k + 3N, l + 2N]*Gl0.val[kPi + 3N, lPj + 2N] - Gl0.val[k + 3N, lPj + 3N]*Gl0.val[kPi + 3N, l + 3N] + 
    Gl0.val[k + 3N, l + 3N]*Gl0.val[kPi + 3N, lPj + 3N])*(-(G0l.val[lPj, k + N]*G0l.val[l + N, kPi]) + (G0l.val[l + N, k + N] - id*I[l, k])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
  (Gl0.val[k, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] + 
    Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*(-(G0l.val[l, k + 2N]*G0l.val[lPj, kPi + N]) + G0l.val[l, kPi + N]*G0l.val[lPj, k + 2N] - 
    G0l.val[l + N, k + 3N]*G0l.val[lPj + N, kPi] + G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 3N] + G0l.val[lPj, k + 3N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + 
    G0l.val[lPj + N, k + 2N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + G0l.val[l, k + 3N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi]) + 
    G0l.val[l + N, k + 2N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
  (Gl0.val[k + N, l + 3N]*Gl0.val[kPi + 2N, lPj + 2N] - Gl0.val[k + N, lPj + 2N]*Gl0.val[kPi + 2N, l + 3N] - Gl0.val[k, l + 3N]*Gl0.val[kPi + 3N, lPj + 2N] + 
    Gl0.val[k, lPj + 2N]*Gl0.val[kPi + 3N, l + 3N])*(-(G0l.val[l, k + 3N]*G0l.val[lPj, kPi + N]) + G0l.val[l, kPi + N]*G0l.val[lPj, k + 3N] - G0l.val[l + N, k + 2N]*G0l.val[lPj + N, kPi] + 
    G0l.val[l + N, kPi]*G0l.val[lPj + N, k + 2N] + G0l.val[lPj, k + 2N]*(-G0l.val[l, kPi] + id*I[l, kPi]) + G0l.val[lPj + N, k + 3N]*(-G0l.val[l + N, kPi + N] + id*I[l, kPi]) + 
    G0l.val[l, k + 2N]*(G0l.val[lPj, kPi] - id*I[lPj, kPi]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) - 
  (Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj] - Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj + N] + Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 2N] - Gl0.val[k, lPj]*Gl0.val[kPi, l + 3N] - 
    Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj] + Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj + N] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 2N] + 
    Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 3N])*(-(G0l.val[l, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) - G0l.val[l + N, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + 
    G0l.val[l + N, kPi + 2N]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k]) + G0l.val[l, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
  (Gl0.val[k, l + 2N]*Gl0.val[kPi, lPj] - Gl0.val[k, l + 3N]*Gl0.val[kPi, lPj + N] - Gl0.val[k, lPj]*Gl0.val[kPi, l + 2N] + Gl0.val[k, lPj + N]*Gl0.val[kPi, l + 3N] - 
    Gl0.val[k + N, l + 2N]*Gl0.val[kPi + N, lPj] + Gl0.val[k + N, l + 3N]*Gl0.val[kPi + N, lPj + N] + Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + 2N] - 
    Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + 3N])*(-(G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 2N, k + 3N]) - G0l.val[l, k + 3N]*G0l.val[lPj + 3N, kPi + 2N] + 
    G0l.val[l, kPi + 2N]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
  (Gl0.val[k, l + N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + N] - Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj] + Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + N])*
   (G0l.val[lPj + 3N, kPi + 2N]*(-G0l.val[l + 3N, k + 3N] + id*I[l, k]) + G0l.val[lPj + 2N, k + 3N]*(-G0l.val[l + 2N, kPi + 2N] + id*I[l, kPi]) + 
    G0l.val[l + 3N, kPi + 2N]*(G0l.val[lPj + 3N, k + 3N] - id*I[lPj, k]) + G0l.val[l + 2N, k + 3N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
  (Gl0.val[kPi, lPj]*Gl0.val[k + N, l] - Gl0.val[kPi, l]*Gl0.val[k + N, lPj] - Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + N] + Gl0.val[kPi, l + N]*Gl0.val[k + N, lPj + N])*
   (G0l.val[lPj + 2N, kPi + 3N]*(-G0l.val[l + 3N, k + 3N] + id*I[l, k]) + G0l.val[lPj + 2N, k + 3N]*(G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi]) + 
    G0l.val[l + 3N, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l + 3N, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
  (-(Gl0.val[k, lPj]*Gl0.val[kPi, l]) + Gl0.val[k, l]*Gl0.val[kPi, lPj] + Gl0.val[k, lPj + N]*Gl0.val[kPi, l + N] - Gl0.val[k, l + N]*Gl0.val[kPi, lPj + N] + Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l] - 
    Gl0.val[k + N, l]*Gl0.val[kPi + N, lPj] - Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + N] + Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj + N])*
   (-(G0l.val[lPj + 2N, k + 3N]*G0l.val[l + 3N, kPi + 2N]) + (G0l.val[l + 3N, k + 3N] - id*I[l, k])*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
  (-(Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj + N] + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 2N] - Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 3N])*
   (G0l.val[l + N, kPi + 3N]*G0l.val[lPj + 2N, k + 3N] - G0l.val[l + N, k + 3N]*G0l.val[lPj + 2N, kPi + 3N] + G0l.val[l, kPi + 2N]*G0l.val[lPj + 3N, k + 2N] - 
    G0l.val[l, k + 2N]*G0l.val[lPj + 3N, kPi + 2N] + G0l.val[l + N, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + 
    G0l.val[l, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + G0l.val[l + N, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi]) + 
    G0l.val[l, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + (-(Gl0.val[kPi, l + 3N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, l + 2N]*Gl0.val[k + N, lPj + N] - 
    Gl0.val[kPi, lPj + N]*Gl0.val[k + N, l + 2N] + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + 3N])*(G0l.val[l, kPi + 3N]*G0l.val[lPj + 2N, k + 3N] - 
    G0l.val[l, k + 3N]*G0l.val[lPj + 2N, kPi + 3N] + G0l.val[l + N, kPi + 2N]*G0l.val[lPj + 3N, k + 2N] - G0l.val[l + N, k + 2N]*G0l.val[lPj + 3N, kPi + 2N] + 
    G0l.val[l, kPi + 2N]*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + G0l.val[l + N, kPi + 3N]*(-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + 
    G0l.val[l, k + 2N]*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi]) + G0l.val[l + N, k + 3N]*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) + 
  (-(Gl0.val[kPi + 2N, l + 3N]*Gl0.val[k + 3N, lPj + 2N]) + Gl0.val[kPi + 2N, lPj + 2N]*Gl0.val[k + 3N, l + 3N])*
   (G0l.val[l, kPi + N]*G0l.val[lPj, k + N] - G0l.val[l, k + N]*G0l.val[lPj, kPi + N] + G0l.val[l + N, kPi]*G0l.val[lPj + N, k] - G0l.val[l + N, k]*G0l.val[lPj + N, kPi] + 
    (G0l.val[l, kPi] - id*I[l, kPi])*(-G0l.val[lPj, k] + id*I[lPj, k]) + (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k]) + 
    (-G0l.val[l, k] + id*I[l, k])*(-G0l.val[lPj, kPi] + id*I[lPj, kPi]) + (-G0l.val[l + N, k + N] + id*I[l, k])*(-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi])) + 
  (-(Gl0.val[kPi, l + N]*Gl0.val[k + N, lPj]) + Gl0.val[kPi, lPj]*Gl0.val[k + N, l + N])*(G0l.val[l + 2N, kPi + 3N]*G0l.val[lPj + 2N, k + 3N] - 
    G0l.val[l + 2N, k + 3N]*G0l.val[lPj + 2N, kPi + 3N] + G0l.val[l + 3N, kPi + 2N]*G0l.val[lPj + 3N, k + 2N] - G0l.val[l + 3N, k + 2N]*G0l.val[lPj + 3N, kPi + 2N] + 
    (G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi])*(-G0l.val[lPj + 2N, k + 2N] + id*I[lPj, k]) + (G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi])*
     (-G0l.val[lPj + 3N, k + 3N] + id*I[lPj, k]) + (-G0l.val[l + 2N, k + 2N] + id*I[l, k])*(-G0l.val[lPj + 2N, kPi + 2N] + id*I[lPj, kPi]) + 
    (-G0l.val[l + 3N, k + 3N] + id*I[l, k])*(-G0l.val[lPj + 3N, kPi + 3N] + id*I[lPj, kPi])));

    tcoscos3=-((Gl0.val[kPi + 2N, l + 2N]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
    Gl0.val[kPi + 3N, l + 2N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]) - 
    Gl0.val[k + 2N, l + 2N]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + Gl0.val[k + 3N, l + 2N]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
    Gl0.val[k, l + 2N]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + Gl0.val[k + N, l + 2N]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
   ((G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj, kPi] + (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*
     G0l.val[l + N, kPi] - (G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi] + (G00.val[l, lPj + N] + G00.val[l + N, lPj])*G0l.val[lPj + 2N, kPi] - 
    (G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 3N, kPi] - (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*
     (G0l.val[l, kPi] - id*I[l, kPi]) - id*(G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*I[lPj, kPi])) - 
 (Gl0.val[kPi + 3N, l + 3N]*(-Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + 2N, l + 3N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]) - 
   Gl0.val[k + 3N, l + 3N]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + Gl0.val[k + 2N, l + 3N]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
   Gl0.val[k + N, l + 3N]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + Gl0.val[k, l + 3N]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
  ((G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l, kPi + N] - 
   (G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi + N] + (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj + N, kPi + N] - 
   (G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 2N, kPi + N] + (G00.val[l, lPj + N] + G00.val[l + N, lPj])*G0l.val[lPj + 3N, kPi + N] + 
   (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) - 
   id*(G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*I[lPj, kPi]) - 
 (Gl0.val[kPi + 2N, l + 3N]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + 3N, l + 3N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]) - 
   Gl0.val[k + 2N, l + 3N]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + Gl0.val[k + 3N, l + 3N]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
   Gl0.val[k, l + 3N]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + Gl0.val[k + N, l + 3N]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
  (-((G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi]) + (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*
    G0l.val[l + N, kPi] + (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj + N, kPi] - (G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 2N, kPi] + 
   (G00.val[l, lPj + N] + G00.val[l + N, lPj])*G0l.val[lPj + 3N, kPi] + (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*
    (G0l.val[l, kPi] - id*I[l, kPi]) + id*(G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*I[lPj, kPi]) - 
 (Gl0.val[kPi + 3N, l + 2N]*(-Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + 2N, l + 2N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]) - 
   Gl0.val[k + 3N, l + 2N]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + Gl0.val[k + 2N, l + 2N]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
   Gl0.val[k + N, l + 2N]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + Gl0.val[k, l + 2N]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
  ((-G00.val[lPj, lPj + 3N] + G00.val[lPj + N, lPj + 2N] - G00.val[lPj + 2N, lPj + N] + G00.val[lPj + 3N, lPj])*G0l.val[l, kPi + N] + 
   (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj, kPi + N] - (G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi + N] + 
   (G00.val[l, lPj + N] + G00.val[l + N, lPj])*G0l.val[lPj + 2N, kPi + N] - (G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 3N, kPi + N] - 
   (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*(G0l.val[l + N, kPi + N] - id*I[l, kPi]) + 
   id*(G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*I[lPj, kPi]) - 
 (-(Gl0.val[k + 2N, l]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k])) + Gl0.val[k + 3N, l]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
   Gl0.val[k, l]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) + Gl0.val[k + N, l]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
   Gl0.val[kPi, l]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + N, l]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
  ((G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*G0l.val[lPj, kPi + 2N] - (G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj + N, kPi + 2N] + 
   (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 2N, kPi + 2N] + (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + 
     G00.val[lPj + 3N, lPj + N])*G0l.val[l + 3N, kPi + 2N] - (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi + 2N] - 
   (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*(G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi]) - 
   id*(G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*I[lPj, kPi]) - 
 (-(Gl0.val[k + 3N, l + N]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k])) + Gl0.val[k + 2N, l + N]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
   Gl0.val[k + N, l + N]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) + Gl0.val[k, l + N]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
   Gl0.val[kPi + N, l + N]*(-Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k]) + 
   Gl0.val[kPi, l + N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]))*
  (-((G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj, kPi + 3N]) + (G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*G0l.val[lPj + N, kPi + 3N] + 
   (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l + 2N, kPi + 3N] - 
   (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi + 3N] + (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 3N, kPi + 3N] + 
   (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*(G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi]) - 
   id*(G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*I[lPj, kPi]) + 
 (-(Gl0.val[k + 3N, l]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k])) + Gl0.val[k + 2N, l]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
   Gl0.val[k + N, l]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) + Gl0.val[k, l]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
   Gl0.val[kPi + N, l]*(-Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k]) + 
   Gl0.val[kPi, l]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]))*
  (-((G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*G0l.val[lPj, kPi + 3N]) + (G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj + N, kPi + 3N] + 
   (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*G0l.val[l + 2N, kPi + 3N] - 
   (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 2N, kPi + 3N] + (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi + 3N] + 
   (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*(G0l.val[l + 3N, kPi + 3N] - id*I[l, kPi]) - 
   id*(G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*I[lPj, kPi]) - 
 (-(Gl0.val[k + 2N, l + N]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k])) + Gl0.val[k + 3N, l + N]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
   Gl0.val[k, l + N]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) + Gl0.val[k + N, l + N]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
   Gl0.val[kPi, l + N]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + N, l + N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
  (-((G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj, kPi + 2N]) + (G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*G0l.val[lPj + N, kPi + 2N] - 
   (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi + 2N] + (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - 
     G00.val[lPj + 3N, lPj])*G0l.val[l + 3N, kPi + 2N] + (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 3N, kPi + 2N] + 
   (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*(G0l.val[l + 2N, kPi + 2N] - id*I[l, kPi]) + 
   id*(G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*I[lPj, kPi]) - 
 (Gl0.val[kPi + 2N, l]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + 3N, l]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]) - Gl0.val[k + 2N, l]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + 
   Gl0.val[k + 3N, l]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - Gl0.val[k, l]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + 
   Gl0.val[k + N, l]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*(-((G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj + N, kPi]) + 
   (-G00.val[lPj, lPj + 3N] + G00.val[lPj + N, lPj + 2N] - G00.val[lPj + 2N, lPj + N] + G00.val[lPj + 3N, lPj])*G0l.val[l + 2N, kPi] + 
   (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 2N, kPi] + (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + 
     G00.val[lPj + 3N, lPj + N])*G0l.val[l + 3N, kPi] - (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi] + 
   (G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
 (Gl0.val[kPi + 2N, l + N]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + 3N, l + N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]) - 
   Gl0.val[k + 2N, l + N]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + Gl0.val[k + 3N, l + N]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
   Gl0.val[k, l + N]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + Gl0.val[k + N, l + N]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
  ((G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*G0l.val[lPj + N, kPi] + (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - 
     G00.val[lPj + 3N, lPj + N])*G0l.val[l + 2N, kPi] - (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi] + 
   (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*G0l.val[l + 3N, kPi] + 
   (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 3N, kPi] - (G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) - 
 (Gl0.val[kPi + 3N, l + N]*(-Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + 2N, l + N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]) - 
   Gl0.val[k + 3N, l + N]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + Gl0.val[k + 2N, l + N]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
   Gl0.val[k + N, l + N]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + Gl0.val[k, l + N]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
  (-((G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*G0l.val[lPj, kPi + N]) + (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - 
     G00.val[lPj + 3N, lPj + N])*G0l.val[l + 2N, kPi + N] - (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 2N, kPi + N] + 
   (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*G0l.val[l + 3N, kPi + N] + 
   (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 3N, kPi + N] + (G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*
    (G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) - (Gl0.val[kPi + 3N, l]*(-Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + 2N, l]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]) - 
   Gl0.val[k + 3N, l]*(Gll.val[kPi + 2N, k + N] + Gll.val[kPi + 3N, k]) + Gl0.val[k + 2N, l]*(Gll.val[kPi + 2N, k] + Gll.val[kPi + 3N, k + N]) - 
   Gl0.val[k + N, l]*(Gll.val[kPi + 2N, k + 3N] + Gll.val[kPi + 3N, k + 2N]) + Gl0.val[k, l]*(Gll.val[kPi + 2N, k + 2N] + Gll.val[kPi + 3N, k + 3N]))*
  ((G00.val[l + 2N, lPj + 3N] + G00.val[l + 3N, lPj + 2N])*G0l.val[lPj, kPi + N] + (-G00.val[lPj, lPj + 3N] + G00.val[lPj + N, lPj + 2N] - G00.val[lPj + 2N, lPj + N] + 
     G00.val[lPj + 3N, lPj])*G0l.val[l + 2N, kPi + N] + (G00.val[l + 2N, lPj + N] + G00.val[l + 3N, lPj])*G0l.val[lPj + 2N, kPi + N] + 
   (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + G00.val[lPj + 3N, lPj + N])*G0l.val[l + 3N, kPi + N] - 
   (G00.val[l + 2N, lPj] + G00.val[l + 3N, lPj + N])*G0l.val[lPj + 3N, kPi + N] - (G00.val[l + 2N, lPj + 2N] + G00.val[l + 3N, lPj + 3N])*
    (G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + (-(Gl0.val[k + 2N, l + 2N]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k])) + 
   Gl0.val[k + 3N, l + 2N]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - Gl0.val[k, l + 2N]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) + 
   Gl0.val[k + N, l + 2N]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + Gl0.val[kPi, l + 2N]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - 
     Gll.val[k + 3N, k]) + Gl0.val[kPi + N, l + 2N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
  ((G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*G0l.val[l, kPi + 2N] - 
   (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj, kPi + 2N] + (G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*
    G0l.val[l + N, kPi + 2N] + (G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi + 2N] + (G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 3N, kPi + 2N] - 
   (G00.val[l, lPj + N] + G00.val[l + N, lPj])*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) - 
 (-(Gl0.val[k + 2N, l + 3N]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k])) + Gl0.val[k + 3N, l + 3N]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
   Gl0.val[k, l + 3N]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) + Gl0.val[k + N, l + 3N]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
   Gl0.val[kPi, l + 3N]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
   Gl0.val[kPi + N, l + 3N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
  ((G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l, kPi + 2N] - 
   (G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi + 2N] + (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*
    G0l.val[l + N, kPi + 2N] + (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj + N, kPi + 2N] + (G00.val[l, lPj + N] + G00.val[l + N, lPj])*G0l.val[lPj + 3N, kPi + 2N] - 
   (G00.val[l, lPj] + G00.val[l + N, lPj + N])*(G0l.val[lPj + 2N, kPi + 2N] - id*I[lPj, kPi])) + 
 (Gl0.val[k + 3N, l + 3N]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k]) - Gl0.val[k + 2N, l + 3N]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) + 
   Gl0.val[k + N, l + 3N]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) - Gl0.val[k, l + 3N]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
   Gl0.val[kPi + N, l + 3N]*(Gll.val[k, k + 3N] - Gll.val[k + N, k + 2N] + Gll.val[k + 2N, k + N] - Gll.val[k + 3N, k]) + 
   Gl0.val[kPi, l + 3N]*(Gll.val[k, k + 2N] - Gll.val[k + N, k + 3N] + Gll.val[k + 2N, k] - Gll.val[k + 3N, k + N]))*
  ((G00.val[lPj, lPj + 2N] - G00.val[lPj + N, lPj + 3N] + G00.val[lPj + 2N, lPj] - G00.val[lPj + 3N, lPj + N])*G0l.val[l, kPi + 3N] - 
   (G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj, kPi + 3N] + (G00.val[lPj, lPj + 3N] - G00.val[lPj + N, lPj + 2N] + G00.val[lPj + 2N, lPj + N] - G00.val[lPj + 3N, lPj])*
    G0l.val[l + N, kPi + 3N] + (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj + N, kPi + 3N] - (G00.val[l, lPj] + G00.val[l + N, lPj + N])*G0l.val[lPj + 2N, kPi + 3N] + 
   (G00.val[l, lPj + N] + G00.val[l + N, lPj])*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi])) - 
 (-(Gl0.val[k + 3N, l + 2N]*(Gll.val[kPi, k + N] + Gll.val[kPi + N, k])) + Gl0.val[k + 2N, l + 2N]*(Gll.val[kPi, k] + Gll.val[kPi + N, k + N]) - 
   Gl0.val[k + N, l + 2N]*(Gll.val[kPi, k + 3N] + Gll.val[kPi + N, k + 2N]) + Gl0.val[k, l + 2N]*(Gll.val[kPi, k + 2N] + Gll.val[kPi + N, k + 3N]) + 
   Gl0.val[kPi + N, l + 2N]*(-Gll.val[k, k + 3N] + Gll.val[k + N, k + 2N] - Gll.val[k + 2N, k + N] + Gll.val[k + 3N, k]) + 
   Gl0.val[kPi, l + 2N]*(-Gll.val[k, k + 2N] + Gll.val[k + N, k + 3N] - Gll.val[k + 2N, k] + Gll.val[k + 3N, k + N]))*
  ((-G00.val[lPj, lPj + 3N] + G00.val[lPj + N, lPj + 2N] - G00.val[lPj + 2N, lPj + N] + G00.val[lPj + 3N, lPj])*G0l.val[l, kPi + 3N] + 
   (G00.val[l, lPj + 3N] + G00.val[l + N, lPj + 2N])*G0l.val[lPj, kPi + 3N] + (-G00.val[lPj, lPj + 2N] + G00.val[lPj + N, lPj + 3N] - G00.val[lPj + 2N, lPj] + 
     G00.val[lPj + 3N, lPj + N])*G0l.val[l + N, kPi + 3N] - (G00.val[l, lPj + 2N] + G00.val[l + N, lPj + 3N])*G0l.val[lPj + N, kPi + 3N] + 
   (G00.val[l, lPj + N] + G00.val[l + N, lPj])*G0l.val[lPj + 2N, kPi + 3N] - (G00.val[l, lPj] + G00.val[l + N, lPj + N])*(G0l.val[lPj + 3N, kPi + 3N] - id*I[lPj, kPi]));

   return tcoscos1 + tcoscos2 -4tcoscos3 
end


@inline Base.@propagate_inbounds function full_B1p_Q1Q2_kernel(
  mc, model::TwoBandModel, klkPlP::NTuple{4}, packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_X_symm
)
  k, l, kPi, lPj = klkPlP   #k, l, k+i, l+j
  G00, G0l, Gl0, Gll = packed_greens
  N = length(lattice(model))
  id = I[G0l.k, G0l.l] 

  tcoscos1=0;
  tcoscos2=4*(conj(Gl0.val[k, lPj + N])*(-(conj(G0l.val[l, kPi + N])*conj(G0l.val[lPj, k + N])) + conj(G0l.val[l, k + N])*conj(G0l.val[lPj, kPi + N]) - 
  G0l.val[l, kPi + N]*G0l.val[lPj, k + N] + G0l.val[l, k + N]*G0l.val[lPj, kPi + N])*Gl0.val[kPi, l + N] + conj(G0l.val[l, kPi + N])*G0l.val[lPj, k + N]*
 (-(conj(Gl0.val[k, lPj + N])*conj(Gl0.val[kPi, l + N])) + conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi, lPj + N]) - Gl0.val[k, lPj + N]*Gl0.val[kPi, l + N] + 
  Gl0.val[k, l + N]*Gl0.val[kPi, lPj + N]) + conj(Gl0.val[kPi + N, l])*(-(conj(G0l.val[l + N, kPi])*conj(G0l.val[lPj + N, k])) + 
  conj(G0l.val[l + N, k])*conj(G0l.val[lPj + N, kPi]) - G0l.val[l + N, kPi]*G0l.val[lPj + N, k] + G0l.val[l + N, k]*G0l.val[lPj + N, kPi])*Gl0.val[k + N, lPj] + 
conj(G0l.val[lPj + N, k])*G0l.val[l + N, kPi]*(-(conj(Gl0.val[k + N, lPj])*conj(Gl0.val[kPi + N, l])) + conj(Gl0.val[k + N, l])*conj(Gl0.val[kPi + N, lPj]) - 
  Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l] + Gl0.val[k + N, l]*Gl0.val[kPi + N, lPj]) - 
(conj(Gl0.val[k + N, l + N])*conj(Gl0.val[kPi + N, lPj]) - conj(Gl0.val[k + N, lPj])*conj(Gl0.val[kPi + N, l + N]) + Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj] - 
  Gl0.val[k + N, lPj]*Gl0.val[kPi + N, l + N])*(G0l.val[lPj + N, kPi]*(conj(G0l.val[l, k]) - id*I[l, k]) + conj(G0l.val[lPj + N, k])*(G0l.val[l, kPi] - id*I[l, kPi])) + 
(-(conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi + N, l])) + conj(Gl0.val[k, l])*conj(Gl0.val[kPi + N, lPj]) - Gl0.val[k, lPj]*Gl0.val[kPi + N, l] + Gl0.val[k, l]*Gl0.val[kPi + N, lPj])*
 (conj(G0l.val[lPj + N, kPi])*(-G0l.val[l + N, k + N] + id*I[l, k]) + G0l.val[l + N, kPi]*(-conj(G0l.val[lPj + N, k + N]) + id*I[lPj, k])) + 
(conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi, l]) - conj(Gl0.val[k, l])*conj(Gl0.val[kPi, lPj]) + Gl0.val[k, lPj]*Gl0.val[kPi, l] - Gl0.val[k, l]*Gl0.val[kPi, lPj])*
 (G0l.val[l + N, k + N] - id*I[l, k])*(-conj(G0l.val[lPj + N, kPi + N]) + id*I[lPj, kPi]) + 
(conj(Gl0.val[k + N, lPj + N])*conj(Gl0.val[kPi + N, l + N]) - conj(Gl0.val[k + N, l + N])*conj(Gl0.val[kPi + N, lPj + N]) + 
  Gl0.val[k + N, lPj + N]*Gl0.val[kPi + N, l + N] - Gl0.val[k + N, l + N]*Gl0.val[kPi + N, lPj + N])*(conj(G0l.val[l, k]) - id*I[l, k])*(-G0l.val[lPj, kPi] + id*I[lPj, kPi]) - 
(conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi, lPj]) - conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi, l + N]) + Gl0.val[k, l + N]*Gl0.val[kPi, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi, l + N])*
 (conj(G0l.val[l, kPi + N])*(G0l.val[lPj + N, k + N] - id*I[lPj, k]) + G0l.val[l, k + N]*(conj(G0l.val[lPj + N, kPi + N]) - id*I[lPj, kPi])) + 
(conj(Gl0.val[k, lPj + N])*conj(Gl0.val[kPi + N, l + N]) - conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi + N, lPj + N]) + Gl0.val[k, lPj + N]*Gl0.val[kPi + N, l + N] - 
  Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj + N])*(G0l.val[lPj, k + N]*(conj(G0l.val[l, kPi]) - id*I[l, kPi]) + conj(G0l.val[l, k + N])*(G0l.val[lPj, kPi] - id*I[lPj, kPi])) + 
conj(Gl0.val[k, l])*Gl0.val[kPi, lPj]*((conj(G0l.val[l + N, kPi + N]) - id*I[l, kPi])*(-conj(G0l.val[lPj + N, k + N]) + id*I[lPj, k]) + 
  (G0l.val[l + N, kPi + N] - id*I[l, kPi])*(-G0l.val[lPj + N, k + N] + id*I[lPj, k]) + (conj(G0l.val[l + N, k + N]) - id*I[l, k])*
   (conj(G0l.val[lPj + N, kPi + N]) - id*I[lPj, kPi]) + (G0l.val[l + N, k + N] - id*I[l, k])*(G0l.val[lPj + N, kPi + N] - id*I[lPj, kPi])) + 
(conj(Gl0.val[k, l + N])*Gl0.val[kPi, lPj] + conj(Gl0.val[k, lPj])*Gl0.val[kPi, l + N])*(conj(G0l.val[l, kPi + N])*conj(G0l.val[lPj + N, k + N]) - 
  conj(G0l.val[l, k + N])*conj(G0l.val[lPj + N, kPi + N]) + G0l.val[l, kPi + N]*G0l.val[lPj + N, k + N] - G0l.val[l, k + N]*G0l.val[lPj + N, kPi + N] + 
  2*id*I[lPj, kPi]*real(G0l.val[l, k + N]) - 2*id*I[lPj, k]*real(G0l.val[l, kPi + N])) + conj(Gl0.val[kPi + N, lPj + N])*Gl0.val[k + N, l + N]*
 (-(conj(G0l.val[l, kPi])*conj(G0l.val[lPj, k])) + conj(G0l.val[l, k])*conj(G0l.val[lPj, kPi]) - G0l.val[l, kPi]*G0l.val[lPj, k] + G0l.val[l, k]*G0l.val[lPj, kPi] + 
  2*id*(-(I[lPj, kPi]*real(G0l.val[l, k])) + I[lPj, k]*real(G0l.val[l, kPi]) + I[l, kPi]*(-(id*I[lPj, k]) + real(G0l.val[lPj, k])) + I[l, k]*(id*I[lPj, kPi] - real(G0l.val[lPj, kPi])))) - 
(conj(Gl0.val[kPi + N, lPj + N])*Gl0.val[k, l + N] + conj(Gl0.val[k, lPj + N])*Gl0.val[kPi + N, l + N])*(conj(G0l.val[l, k + N])*conj(G0l.val[lPj, kPi]) - 
  conj(G0l.val[l, kPi])*conj(G0l.val[lPj, k + N]) + G0l.val[l, k + N]*G0l.val[lPj, kPi] - G0l.val[l, kPi]*G0l.val[lPj, k + N] - 2*id*I[lPj, kPi]*real(G0l.val[l, k + N]) + 
  2*id*I[l, kPi]*real(G0l.val[lPj, k + N])) - (conj(Gl0.val[kPi + N, l])*Gl0.val[k, lPj] + conj(Gl0.val[k, l])*Gl0.val[kPi + N, lPj])*
 (conj(G0l.val[l + N, k + N])*conj(G0l.val[lPj + N, kPi]) - conj(G0l.val[l + N, kPi])*conj(G0l.val[lPj + N, k + N]) + G0l.val[l + N, k + N]*G0l.val[lPj + N, kPi] - 
  G0l.val[l + N, kPi]*G0l.val[lPj + N, k + N] + 2*id*I[lPj, k]*real(G0l.val[l + N, kPi]) - 2*id*I[l, k]*real(G0l.val[lPj + N, kPi])) + 
(conj(Gl0.val[kPi + N, l + N])*Gl0.val[k + N, lPj] + conj(Gl0.val[kPi + N, lPj])*Gl0.val[k + N, l + N])*(conj(G0l.val[l, kPi])*conj(G0l.val[lPj + N, k]) - 
  conj(G0l.val[l, k])*conj(G0l.val[lPj + N, kPi]) + G0l.val[l, kPi]*G0l.val[lPj + N, k] - G0l.val[l, k]*G0l.val[lPj + N, kPi] - 2*id*I[l, kPi]*real(G0l.val[lPj + N, k]) + 
  2*id*I[l, k]*real(G0l.val[lPj + N, kPi])) + (conj(Gl0.val[k, l + N])*conj(Gl0.val[kPi + N, lPj]) - conj(Gl0.val[k, lPj])*conj(Gl0.val[kPi + N, l + N]) + 
  Gl0.val[k, l + N]*Gl0.val[kPi + N, lPj] - Gl0.val[k, lPj]*Gl0.val[kPi + N, l + N])*(conj(G0l.val[lPj + N, k + N])*G0l.val[l, kPi] + conj(G0l.val[lPj + N, kPi])*G0l.val[l, k + N] + 
  conj(G0l.val[l, k + N])*G0l.val[lPj + N, kPi] + conj(G0l.val[l, kPi])*G0l.val[lPj + N, k + N] - 2*id*I[lPj, k]*real(G0l.val[l, kPi]) + 
  2*id*I[l, kPi]*(id*I[lPj, k] - real(G0l.val[lPj + N, k + N]))) + (conj(Gl0.val[kPi + N, l + N])*Gl0.val[k, lPj] + conj(Gl0.val[kPi + N, lPj])*Gl0.val[k, l + N] + 
  conj(Gl0.val[k, l + N])*Gl0.val[kPi + N, lPj] + conj(Gl0.val[k, lPj])*Gl0.val[kPi + N, l + N])*(conj(G0l.val[l, k + N])*conj(G0l.val[lPj + N, kPi]) - 
  conj(G0l.val[l, kPi])*conj(G0l.val[lPj + N, k + N]) + G0l.val[l, k + N]*G0l.val[lPj + N, kPi] - G0l.val[l, kPi]*G0l.val[lPj + N, k + N] + 
  2*id*(I[lPj, k]*real(G0l.val[l, kPi]) + I[l, kPi]*(-(id*I[lPj, k]) + real(G0l.val[lPj + N, k + N])))));

  tcoscos3=2*(4*imag(G00.val[lPj, lPj + N] - G00.val[lPj + N, lPj])*
  imag((G0l.val[l + N, kPi + N] - id*I[l, kPi])*(conj(Gl0.val[k + N, l])*real(Gll.val[kPi, k]) - conj(Gl0.val[k, l])*real(Gll.val[kPi, k + N])) + 
    G0l.val[l, kPi + N]*(-(conj(Gl0.val[k + N, l + N])*real(Gll.val[kPi, k])) + conj(Gl0.val[k, l + N])*real(Gll.val[kPi, k + N])) + 
    (G0l.val[l, kPi] - id*I[l, kPi])*(conj(Gl0.val[k + N, l + N])*real(Gll.val[kPi + N, k]) - conj(Gl0.val[k, l + N])*real(Gll.val[kPi + N, k + N])) + 
    conj(G0l.val[l + N, kPi])*(Gl0.val[k + N, l]*real(Gll.val[kPi + N, k]) - Gl0.val[k, l]*real(Gll.val[kPi + N, k + N]))) + 
 2*imag(Gll.val[k, k + N] - Gll.val[k + N, k])*
  (imag(2*Gl0.val[kPi, l + N]*((-conj(G0l.val[lPj + N, kPi + N]) + id*I[lPj, kPi])*real(G00.val[l, lPj]) + conj(G0l.val[lPj, kPi + N])*real(G00.val[l, lPj + N])) + 
     2*Gl0.val[kPi + N, l + N]*(conj(G0l.val[lPj + N, kPi])*real(G00.val[l, lPj]) + (-conj(G0l.val[lPj, kPi]) + id*I[lPj, kPi])*real(G00.val[l, lPj + N])) - 
     Gl0.val[kPi, l]*(-2*(conj(G0l.val[lPj + N, kPi + N]) - id*I[lPj, kPi])*real(G00.val[l + N, lPj]) + 2*conj(G0l.val[lPj, kPi + N])*real(G00.val[l + N, lPj + N])) - 
     2*Gl0.val[kPi + N, l]*(conj(G0l.val[lPj + N, kPi])*real(G00.val[l + N, lPj]) + (-conj(G0l.val[lPj, kPi]) + id*I[lPj, kPi])*real(G00.val[l + N, lPj + N]))) + 
   2*imag(G00.val[lPj, lPj + N] - G00.val[lPj + N, lPj])*real(-(conj(Gl0.val[kPi + N, l + N])*G0l.val[l, kPi]) + conj(Gl0.val[kPi, l + N])*G0l.val[l, kPi + N] + 
      conj(Gl0.val[kPi + N, l])*G0l.val[l + N, kPi] - conj(Gl0.val[kPi, l])*G0l.val[l + N, kPi + N] + id*(Gl0.val[kPi, l] + Gl0.val[kPi + N, l + N])*I[l, kPi])) + 
 real((-2*G0l.val[lPj + N, kPi]*real(G00.val[l, lPj]) + 2*(G0l.val[lPj, kPi] - id*I[lPj, kPi])*real(G00.val[l, lPj + N]))*(2*conj(Gl0.val[k + N, l + N])*real(Gll.val[kPi + N, k]) - 
     2*conj(Gl0.val[k, l + N])*real(Gll.val[kPi + N, k + N])) + 2*(((-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi])*real(G00.val[l + N, lPj]) + G0l.val[lPj, kPi + N]*real(G00.val[l + N, lPj + N]))*
      (2*conj(Gl0.val[k + N, l])*real(Gll.val[kPi, k]) - 2*conj(Gl0.val[k, l])*real(Gll.val[kPi, k + N])) - 
     ((-G0l.val[lPj + N, kPi + N] + id*I[lPj, kPi])*real(G00.val[l, lPj]) + G0l.val[lPj, kPi + N]*real(G00.val[l, lPj + N]))*(2*conj(Gl0.val[k + N, l + N])*real(Gll.val[kPi, k]) - 
       2*conj(Gl0.val[k, l + N])*real(Gll.val[kPi, k + N])) + (G0l.val[lPj + N, kPi]*real(G00.val[l + N, lPj]) + (-G0l.val[lPj, kPi] + id*I[lPj, kPi])*real(G00.val[l + N, lPj + N]))*
      (2*conj(Gl0.val[k + N, l])*real(Gll.val[kPi + N, k]) - 2*conj(Gl0.val[k, l])*real(Gll.val[kPi + N, k + N])))));
      
  return tcoscos1 + tcoscos2 -4tcoscos3 
end






