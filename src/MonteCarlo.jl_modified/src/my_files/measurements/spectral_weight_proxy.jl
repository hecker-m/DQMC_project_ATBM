
function spectral_weight_proxy(
    mc::DQMC, model::Model;  greens_iterator=TimeIntegral(mc),
        lattice_iterator = EachSitePairByDistance(), 
        flavor_iterator = FlavorIterator(mc, 1),
        #returns constant number of fermion flavors, i.e. 4
        kernel = my_spec_weight_kernel,
        kwargs...
    )
    return DQMCMeasurement(mc, model, greens_iterator, lattice_iterator,
     flavor_iterator, kernel; kwargs...)
end


#=
@inline Base.@propagate_inbounds function my_spec_weight_kernel(
    mc, model, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv)
    nslices=mc.parameters.slices
    G00, G0l, Gl0, Gll = packed_greens
    if Gll.l==div(nslices,2)
        i, j = ij
        N = length(lattice(mc))
        return Gl0.val[i+(flv-1)*N,j+(flv-1)*N]
   end
end
=#
@inline Base.@propagate_inbounds function my_spec_weight_kernel(
    mc, model, ij::NTuple{2}, Gl0::GreensMatrix, flv)
        i, j = ij
        N = length(lattice(mc))
    return Gl0.val[i+(flv-1)*N,j+(flv-1)*N]  
end

function apply!(
    temp::Array, iter::EachSitePairByDistance, 
    measurement::DQMCMeasurement{typeof(my_spec_weight_kernel)}, mc::DQMC, 
    packed_greens, weight = 1.0)
    G00, G0l, Gl0, Gll = packed_greens
    nslices=mc.parameters.slices
    if Gll.l==div(nslices,2)
        l = lattice(mc)
        Bsrctrg2dir = l[:Bravais_srctrg2dir]::Matrix{Int}
        B = length(unitcell(l)) #in TwoBandModel B=1
        N = length(Bravais(l))  #in TwoBandModel N=L^2

        @inbounds @fastmath for σ in measurement.flavor_iterator    #for flavor iterator being a constant, there is no iteration over it
            for b2 in 1:B, b1 in 1:B
                uc1 = N * (b1-1)
                uc2 = N * (b2-1)
                for trg in 1:N
                    @simd for src in 1:N
                        dir = Bsrctrg2dir[src, trg]
                        temp[dir, b1, b2] += 1.0 * 
                        measurement.kernel(mc, mc.model, (src + uc1, trg + uc2), Gl0, σ)
                    end
                end
            end
        end

    end
    return
end