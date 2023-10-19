
include("/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/MonteCarlo.jl_modified/src/my_files/q_discretization.jl")


function superfluid_stiffness(
    mc::DQMC,  ccs::DQMCMeasurement; 
    shift_dir = [1., 0.])
    β=mc.parameters.beta
    normalize!(shift_dir)
    l=lattice(mc)
    # Find index corresponding to the smallest jump in q-space parallel to shift_dir
    idx_longi = let
        qs = cached_reciprocal_discretization(l)::Matrix{Vector{Float64}}
        ortho = [-shift_dir[2], shift_dir[1]]
        # idx2dist = map(i -> i => dot(qs[i], ortho), eachindex(qs))
        # filter!(t -> abs(t[2]) < 1e-6, idx2dist)
        # minimum(last, idx2dist)[1]

        idx2dist = map(i -> i => [dot(qs[i], ortho), norm(qs[i])], eachindex(qs))
        filter!(t -> abs(t[2][1]) < 1e-6, idx2dist)
        sort!(idx2dist, by= x -> x[2][2])
        idx2dist[2][1]
    end
    idx_trans = let
        qs = cached_reciprocal_discretization(l)::Matrix{Vector{Float64}}
        idx2dist = map(i -> i => [dot(qs[i], shift_dir), norm(qs[i])], eachindex(qs))
        filter!(t -> abs(t[2][1]) < 1e-6, idx2dist)
        sort!(idx2dist, by= x -> x[2][2])
        idx2dist[2][1]
    end
    Λxx = cached_para_ccc(mc, ccs, shift_dir)
    #print(Λxx)
    @show [idx_longi, idx_trans]
    return [0.25 * (Λxx[idx_longi] - Λxx[idx_trans]), β*π/2 *0.25 * (Λxx[idx_longi] - Λxx[idx_trans])], 
        [Λxx[idx_longi]], Λxx[idx_trans], Λxx[1]
end


function superfluid_stiffness(
    mc::DQMC, G::DQMCMeasurement, ccs::DQMCMeasurement; 
    shift_dir = [1., 0.])
    β=mc.parameters.beta
    normalize!(shift_dir)
    l=lattice(mc)
    # Find index corresponding to the smallest jump in q-space parallel to shift_dir
    idx_longi = let
        qs = cached_reciprocal_discretization(l)::Matrix{Vector{Float64}}
        ortho = [-shift_dir[2], shift_dir[1]]
        # idx2dist = map(i -> i => dot(qs[i], ortho), eachindex(qs))
        # filter!(t -> abs(t[2]) < 1e-6, idx2dist)
        # minimum(last, idx2dist)[1]

        idx2dist = map(i -> i => [dot(qs[i], ortho), norm(qs[i])], eachindex(qs))
        filter!(t -> abs(t[2][1]) < 1e-6, idx2dist)
        sort!(idx2dist, by= x -> x[2][2])
        idx2dist[2][1]
    end
    idx_trans = let
        qs = cached_reciprocal_discretization(l)::Matrix{Vector{Float64}}
        idx2dist = map(i -> i => [dot(qs[i], shift_dir), norm(qs[i])], eachindex(qs))
        filter!(t -> abs(t[2][1]) < 1e-6, idx2dist)
        sort!(idx2dist, by= x -> x[2][2])
        idx2dist[2][1]
    end
    Kx  = dia_K(mc, G, shift_dir)
    Λxx = cached_para_ccc(mc, ccs, shift_dir)
    #print(Λxx)
    @show [idx_longi, idx_trans]
    return [0.25 * (Λxx[idx_longi] - Λxx[idx_trans]), β*π/2 *0.25 * (Λxx[idx_longi] - Λxx[idx_trans])], 
        [Λxx[idx_longi], Kx], Λxx[idx_trans], Λxx[1]
end


"""
dia_K(mc, key::Symbol, dir::Vector)
dia_K(mc, m::DQMCMeasurement, dir::Vector)
dia_K(mc, G::Matrix, dir::Vector)

Computes the diamangetic contribution of electromagnetic response of the system
along a given direction `dir`. 
"""
dia_K(mc::DQMC, key::Symbol, shift_dir) = dia_K(mc, mean(mc[key]), shift_dir)
dia_K(mc::DQMC, m::DQMCMeasurement, shift_dir) = dia_K(mc, mean(m), shift_dir)
function dia_K(mc::DQMC, G::AbstractMatrix, shift_dir)

    # Get hoppings
    if !isdefined(mc.stack, :hopping_matrix)
        MonteCarlo.init_hopping_matrices(mc, mc.model)
    end
    T = Matrix(mc.stack.hopping_matrix)

    # Filter directions appropriate to shift_dir
    normalize!(shift_dir) # to be sure
    dirs = directions(lattice(mc))
    idxs = filter(MonteCarlo.hopping_directions(lattice(mc))) do i
        dot(dirs[i], shift_dir) > 1e-6
    end

    # Lattice iterator to map directional indices to (src, trg) pairs
    dir2srctrg = lattice(mc)[:dir2srctrg]
    N = length(lattice(mc))

    Kx = 0.0
        # sum over valid hopping directions
        for i in idxs

            weight = dot(dirs[i], shift_dir)^2

            for (src, trg) in dir2srctrg[i]     #(src, trg) =(r, r+b)
                # c_j^† c_i = δ_ij - G[i, j], but δ is always 0 because this 
                # excludes on-site. See (2)
                # Reverse directions are filtered out before this and explicitly
                # included again here.

                Kx -= weight * kx_kernel_post(mc, G, T, (src, trg), MonteCarlo.field(mc))
            end
        end
    return Kx / N
end
function kx_kernel_post(mc, G, T, st::NTuple{2}, ::AbstractMagnBosonField)
    N = length(lattice(mc))
    r, rPb =st
    return G[rPb, r]*T[r, rPb] + G[N + rPb, N + r]*T[N + r, N + rPb] + 
        G[2N + rPb, 2N + r]*T[2N + r, 2N + rPb] + G[3N + rPb, 3N + r]*T[3N + r, 3N + rPb] + 
        G[r, rPb]*T[rPb, r] + G[N + r, N + rPb]*T[N + rPb, N + r] + 
        G[2N + r, 2N + rPb]*T[2N + rPb, 2N + r] + G[3N + r, 3N + rPb]*T[3N + rPb, 3N + r]
end
function kx_kernel_post(mc, G, T, st::NTuple{2}, ::Union{Discrete_MBF1_X_symm, Discrete_MBF2_symm})
    N = length(lattice(mc))
    r, rPb =st
    return 2*real(G[rPb, r]*T[r, rPb] + G[N + rPb, N + r]*T[N + r, N + rPb] + 
        G[r, rPb]*T[rPb, r] + G[N + r, N + rPb]*T[N + rPb, N + r])
end



# This maybe usefull to deal with different `iter.directions`
_mapping(dirs::Union{Vector, Tuple}) = copy(dirs)

"""
cached_para_ccc(mc, key::Symbol, dir)
cached_para_ccc(mc, m::DQMCMeasurement, dir)
cached_para_ccc(lattice, iter::AbstractLatticeIterator, ccs, dir)

Returns the Fourier transformed current current correlation along a given 
direction. This is the paramagnetic contribution of the electromagnetic response.
"""
function cached_para_ccc(mc::DQMC, key::Symbol, shift_dir; kwargs...)
return cached_para_ccc(mc, mc[key], shift_dir; kwargs...)
end
function cached_para_ccc(mc::DQMC, m::DQMCMeasurement, shift_dir; kwargs...)
return cached_para_ccc(lattice(mc), m.lattice_iterator, mean(m), shift_dir; kwargs...)
end


function reverse_bond_table(l)
# This generates a list reverse_bond[bond_idx]
bs = l.unitcell.bonds
table = zeros(Int, length(bs))

for i in eachindex(table)
    if table[i] == 0
        idx = findfirst(bs) do b
            b.from == bs[i].to &&
            b.to   == bs[i].from &&
            b.uc_shift == .- bs[i].uc_shift
        end::Int

        table[i] = idx
        table[idx] = i
    end
end

return table
end

function cached_para_ccc(l::Lattice, iter::EachBondPairByBravaisDistance, ccs::Array, shift_dir; skip_check = false)
    # Equivalent to the above with a more straight forward lattice iterator

    # Get caches
    equivalency = get!(l, :reverse_bond_table, reverse_bond_table)
    fft_cache = if haskey(l.cache.cache, :fft_cache)
        l.cache.cache[:fft_cache]
    else
        l.cache.cache[:fft_cache] = Matrix{ComplexF64}(undef, l.Ls)
    end::Matrix{ComplexF64}
    qs = cached_reciprocal_discretization(l)::Matrix{Vector{Float64}}

    uc = l.unitcell
    bs = uc.bonds
    v1, v2 = MonteCarlo.lattice_vectors(l)

    # directions of all bonds
    ds = map(bs) do b
        uc.sites[b.to] - uc.sites[b.from] + v1 * b.uc_shift[1] + v2 * b.uc_shift[2]
    end

    # Bond centers relative to the Bravais lattice position
    # !!! Note: The first `+` below could potentially be a minus.
    cs = map(bs) do b
        0.5 * (uc.sites[b.to] + uc.sites[b.from] + v1 * b.uc_shift[1] + v2 * b.uc_shift[2])
    end

    # applicable indices:
    # 1. index into measurement Array (ccs)
    # 2. index into distances ds
    # 3. index into bond centers cs
    applicable = let
        applicable = filter(i -> dot(ds[i], shift_dir) > 0, eachindex(ds))
        if iter.bond_idxs == Colon()
            # All bonds are included, so the mapping is just i -> i
            tuple.(applicable, applicable, applicable)
        else
            # A subset of bonds is included. Typically this means only one of 
            # i -> j and j -> i is considered.
            # We need to find each bond with a component in shift_dir in the 
            # measured subset. Let's say we need the bond i -> j and find
            # 1. the matching bond i -> j. In this case we pass its index in 
            #    bond_idxs (in the measurement) and the index in ds which 
            #    matches the index in cs.
            # 2. the reverse bond j -> i. Here we can recover the measurement 
            #    for i -> j with a factor -1, which we get by considering the 
            #    direction of j -> i instead of i -> j. So we return the index 
            #    into bond_idxs of the reverse bond, the index of j -> i in bs 
            #    for directions ds and the index of i -> j in bs for centers cs.
            # Remark: In fact, cs also has to be inverted in the 2nd case!
            output = NTuple{3, Int}[]

            for idx in applicable
                i = findfirst(isequal(idx), iter.bond_idxs)
                if i === nothing
                    i = findfirst(isequal(equivalency[idx]), iter.bond_idxs)
                    if i === nothing
                        if !skip_check
                            # TODO would be better to also check if hopping 
                            # matrix is 0 for this bond. It should generally 
                            # be ignored then.
                            error("Missing bond $idx (or $(equivalency[idx]))")
                        else
                            continue
                        end
                    end
                    push!(output, (i, equivalency[idx], equivalency[idx]))
                else
                    # push!(output, i => idx)
                    push!(output, (i, idx, idx))
                end
            end

            output
        end
    end
    @show applicable
    # Output
    Λ = zeros(ComplexF64, size(ccs, 1), size(ccs, 2))

    # for (i, bi) in applicable, (j, bj) in applicable
    for (i, di, ci) in applicable, (j, dj, cj) in applicable
        # Fourier transform on Bravais lattice for on pair of bonds (i, j)
        fft_cache[:, :] .= ccs[:, :, i, j]
        fft!(fft_cache)
        
        # weights from bond centers and directions
        weight =  dot(ds[di], shift_dir) * dot(ds[dj], shift_dir)
        for (idx, q) in enumerate(qs)
            # Λ[idx] += fft_cache[idx] * weight * cis(-dot(cs[bi] - cs[bj], q))
            q=q; #no idea why I need this now for longitudinal
            if idx==2
                @show [q, fft_cache[idx], weight, dot(cs[ci] - cs[cj], q)]
            end
            Λ[idx] += fft_cache[idx] * weight * cis(-dot(cs[ci] - cs[cj], q))
        end
    end

    return Λ
end

