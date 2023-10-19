
# TODO generalize
function reciprocal_discretization(l::Lattice{2}, center = true)
    Ls = size(l)
    vecs = MonteCarlo.reciprocal_vectors(l)
    qs = [zeros(2) for _ in 1:Ls[1], _ in 1:Ls[2]]
    
    if center
        temp = zeros(length(Ls))
        for i in 1:Ls[1], j in 1:Ls[2]

            min = Inf
            for n in -1:1, m in -1:1
                temp .= (n + (i-1)/Ls[1]) * vecs[1] .+ (m + (j-1)/Ls[2]) * vecs[2]

                if norm(temp) < min
                    min = norm(temp)
                    qs[i, j] .= temp
                end
            end
        end
    else
        for i in 1:Ls[1], j in 1:Ls[2]
            qs[i, j] = (i-1) / Ls[1] * vecs[1] .+ (j-1) / Ls[2] * vecs[2]
        end
    end

    return qs
end


"""
    fourier(mc, key[, weights])

Needs verification
"""
function fourier(mc::DQMC, key::Symbol, args...)
    measurement = mc[key]
    data = mean(measurement)
    l = lattice(mc)
    qs = reciprocal_distretization(l)
    dirs = directions(l, measurement.lattice_iterator)

    fourier(qs, dirs, data, args...)
end

# TODO
# - The boundschecks below maybe broken, didn't check or think about it yet
# - versions that only transform Bravais lattice vectors/directions (default?)

function fourier(qs::Vector, dirs::Vector, values::Vector)
    @boundscheck length(dirs) == length(qs) == length(values)
    map(qs) do q
        sum(cis(dot(q, v)) * o for (v, o) in zip(dirs, values))
    end
end

# [basis, basis, dir] indexing
function fourier(qs::Vector, dirs::Array{<: Vector, 3}, values::Array{<: Number, 3})
    @boundscheck length(dirs) == length(qs) == length(values)
    map(qs) do q
        sum(cis(dot(q, v)) * o for (v, o) in zip(dirs, values))
    end
end

# [basis, basis, dir, subdir1, subdir2]
function fourier(qs::Vector, dirs::Tuple{Array{<: Vector, 3}, <: Vector}, values::Array{<: Number, 5})
    @boundscheck size(dirs[1]) == size(values)[1:3]
    output = Array{ComplexF64, 3}(undef, length(qs), size(values, 4), size(values, 5))
    
    # @inbounds
    for k in axes(values, 5)
        for j in axes(values, 4)
            vals = view(values, :, :, :, j, k)
            for (i, q) in enumerate(qs)
                output[i, j, k] = sum(cis(dot(q, v)) * o for (v, o) in zip(dirs[1], vals))
            end
        end
    end

    return output
end

# [basis, basis, dir, synced_dir]
function fourier(qs::Vector, dirs::Tuple{Array{<: Vector, 3}, <: Vector}, values::Array{<: Number, 4})
    @boundscheck size(dirs[1]) == size(values)[1:3]
    output = Array{ComplexF64, 2}(undef, length(qs), size(values, 4))

    # @inbounds
    for j in axes(values, 4)
        vals = view(values, :, :, :, j)
        for (i, q) in enumerate(qs)
            output[i, j] = sum(cis(dot(q, v)) * o for (v, o) in zip(dirs[1], vals))
        end
    end

    return output
end

################################################################################
### Cached
################################################################################


function Bravais_dir2indices(lattice::Lattice{2}, ϵ = 1e-6)
    wrap = generate_combinations(Bravais(lattice))
    directions = Vector{Float64}[]
    idxs = Tuple{Int, Int}[]
    sizehint!(directions, length(Bravais(lattice)))

    d = zeros(2)
    new_d = zeros(2)
    p0 = zeros(2)
    p = zeros(2)
    v1, v2 = lattice_vectors(lattice)
    
    for i in 0:lattice.Ls[1]-1
        for j in 0:lattice.Ls[2]-1
            @. p = i * v1 + j * v2
            _apply_wrap!(d, p, p0, wrap, new_d, ϵ)

            # search for d in directions
            # if not present push it, otherwise continue with next iteration
            for (idx, dir) in enumerate(directions)
                new_d .= dir .- d
                b = true
                for v in new_d
                    b = b && (abs(v) < ϵ)
                end
                if b # all_below(new_d, ϵ)
                    error("Every direction must have a unique target site.")
                    # push!(idxs[idx], (i, j))
                    # @goto loop_end
                end
            end
            push!(directions, copy(d))
            push!(idxs, (i+1, j+1))

            # @label loop_end
        end
    end
    
    order = sortperm(directions, by = v -> directed_norm2(v, ϵ))
    return idxs[order]
end

cached_directions(l::Lattice) = get!(l, :dirs, directions)
cached_directions(l::Bravais) = get!(l.l, :dirs, l -> directions(Bravais(l)))

@deprecate generate_reciprocal_discretization(l) reciprocal_discretization(l, false) false

function cached_reciprocal_discretization(l::Lattice{2})
    # The should effectively be the same as their difference is the Periodicity
    # of the reciprocal lattice. 
    # get!(l, :reciprocal_discretization, generate_reciprocal_discretization)
    get!(l, :reciprocal_discretization, reciprocal_discretization)
end