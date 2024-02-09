
###########################
### full Greens function measurement
###########################
"""
    total_greens_measurement(mc::DQMC, model; kwargs...)

    Measures the total Green's function over all time slices, stored in a large Matrix\n
    {τ ∈[0,`Nslices`+1], ℓ∈[1,4] alas (G00, G0l, Gl0, Gll), i ∈[1, `N`*`flv`], j∈[1, `N`*`flv`]}
    \
"""
function total_greens_measurement(
        mc::DQMC, model::Model;
        greens_iterator= TimeIntegral(mc),
        lattice_iterator = total_greens_li(mc), wrapper = nothing, 
        flavor_iterator = nothing,
        kernel = heat_cap_h2_kernel,
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(mc, model, greens_iterator, li, flavor_iterator, kernel; kwargs...)
end

################################################################################
### total greens function lattice iterator
################################################################################
struct total_greens_li <: DirectLatticeIterator
    greens_size::Int   
    N_slices::Int   
end

function total_greens_li(mc::DQMC)
    N = length(lattice(mc))
    flv = unique_flavors(mc)
    Nτ = nslices(mc)
    return total_greens_li(N*flv, Nτ)
end

function output_size(iter::total_greens_li, l::Lattice)
    return (iter.N_slices +1, 4, iter.greens_size , iter.greens_size)
end

function apply!( temp::Array, iter::total_greens_li, 
        measurement, mc::DQMC, packed_greens, weight = 1.0 )

    G00, G0l, Gl0, Gll = packed_greens
    τ = G0l.l 

    temp[τ+1, 1, :, :] .= G00.val[:,:]    
    temp[τ+1, 2, :, :] .= G0l.val[:,:]    
    temp[τ+1, 3, :, :] .= Gl0.val[:,:]    
    temp[τ+1, 4, :, :] .= Gll.val[:,:]    

    return 
end

@inline function finalize_temp!(::total_greens_li, m, mc)
    return nothing    
end
@inline function commit!(::total_greens_li, m) 
    push!(m.observable, m.temp)
end

