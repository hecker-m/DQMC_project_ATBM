"""
    occupation(mc, model, dir; kwargs...)

Generates a measurements of the occupation number.

## Optional Keyword Arguments

- `kernel = occupation_kernel` sets the function representing the Wicks expanded 
expectation value of the measurement. See `occupation_kernel` for details.
- `lattice_iterator = nothing` handled internally
- `flavor_iterator = nothing` handled internally
- kwargs from `DQMCMeasurement`
"""
function occupation(
        mc::DQMC, model::Model; 
        greens_iterator = Greens(),
        lattice_iterator = nothing,
        flavor_iterator = nothing,
        kernel = occupation_kernel,
        capacity = _default_capacity(mc), eltype = geltype(mc),
        obs = begin 
            N = length(lattice(model)) * unique_flavors(mc);
            LogBinner(zeros(geltype(mc), N), capacity=capacity)
            end,
        kwargs...
    )
    eltype = geltype(mc)
    N = length(lattice(model)) * unique_flavors(mc)
    #obs = FullBinner(Vector{eltype})
    temp = Vector{eltype}(undef, N)
    return DQMCMeasurement(
        mc, model, greens_iterator, lattice_iterator, flavor_iterator, kernel; 
        obs = obs, temp = temp, kwargs...
    )
end
"""
    occupation_summed(mc, model, dir; kwargs...)

Same as `occupation(mc, model, dir; kwargs...)` but returns directly a scalar, to save memory in a `FullBinner`. 
"""
function occupation_summed(
    mc::DQMC, model::Model; 
    greens_iterator = Greens(),
    lattice_iterator = nothing,
    flavor_iterator = nothing,
    kernel = occupation_kernel_summed,
    eltype = geltype(mc),
    obs = FullBinner(eltype),
    temp = Vector{geltype(mc)}(undef, 1),
    kwargs...
)
return DQMCMeasurement(
    mc, model, greens_iterator, lattice_iterator, flavor_iterator, kernel; 
    obs = obs, temp = temp, kwargs...
)
end
# To make this work with good efficency we need a custom `measure!` because all 
# the flavor iterators imply sums...

"""
    occupation_kernel(site_index, flavor_index, N_sites, greens_matrix)

Returns `1 - G[i, i] = ⟨nᵢ⟩`.

Note that is expected to be used with `occupation` which generate its own 
lattice and flavor loops internally. 
"""
@inline Base.@propagate_inbounds function occupation_kernel(i, flv, N, G::_GM{T}, mc) where T
    return occupation_kernel(i, flv, N, G, field(mc))
end
@inline Base.@propagate_inbounds function occupation_kernel(i, flv, N, G::_GM{<: Matrix}, ::AbstractField)
    shift = N * (flv - 1)
    return 1 - G.val[i+shift, i+shift]
end


@inline Base.@propagate_inbounds function occupation_kernel(i, flv, N, G::_GM{<: DiagonallyRepeatingMatrix}, ::AbstractField)
    return 1 - G.val.val[i, i]
end

@inline Base.@propagate_inbounds function occupation_kernel(i, flv, N, G::_GM{<: BlockDiagonal}, ::AbstractField)
    return 1 - G.val.blocks[flv][i, i]
end

@inline Base.@propagate_inbounds function occupation_kernel_summed(i, flv, N, G::_GM{T}, mc) where T
    return occupation_kernel(i, flv, N, G, field(mc))
end

@bm function measure!(::Nothing, m::DQMCMeasurement{typeof(occupation_kernel)}, mc::DQMC, packed_greens)
    i = 1
    N = length(lattice(mc))
    @inbounds @fastmath for flv in 1:unique_flavors(mc)
        for n in eachindex(lattice(mc))
            m.temp[i] = m.kernel(n, flv, N, packed_greens, mc)
            i += 1
        end
    end
    push!(m.observable, m.temp)
    nothing
end

@bm function measure!(::Nothing, m::DQMCMeasurement{typeof(occupation_kernel_summed)}, mc::DQMC, packed_greens)
    N = length(lattice(mc))
    m.temp .= zero(geltype(mc))
    @inbounds @fastmath for flv in 1:unique_flavors(mc)
        for i in eachindex(lattice(mc))
            m.temp[1] += m.kernel(i, flv, N, packed_greens, mc)
        end
    end
    push!(m.observable, real(m.temp[1])/N)
    nothing
end