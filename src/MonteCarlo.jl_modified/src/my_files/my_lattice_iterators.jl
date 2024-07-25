
################################################################################
### iterator for a time-averaged energy measurement --- higher accuracy 
### necessary for heat capacity measurement
################################################################################

struct energy_iterator <: DirectLatticeIterator
end

function output_size(iter::energy_iterator, l::Lattice)
    return (1,)
end

function apply!(
    temp::Array, iter::energy_iterator,
    measurement, mc::DQMC, packed_greens, weight=1.0)
    @timeit_debug "apply!(::energy_iterator, ::$(typeof(measurement.kernel)))" begin

        temp[1] += weight * measurement.kernel(mc, mc.model, nothing, packed_greens, nothing)
    end
    return
end

@inline function finalize_temp!(::energy_iterator, m, mc)
    β = mc.parameters.beta
    m.temp .*= (1 / (length(lattice(mc)) * β))  # * 1/N *1/β
end
@inline function commit!(::energy_iterator, m)
    push!(m.observable, m.temp[1])
end




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
output_size(::EachSitePair_B1, l::Lattice) = (1,)

_length(::EachSitePair_B1, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_B1, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachSitePair_B1, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        srcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L = lat.Ls[1]
        @inbounds @fastmath for i in eachindex(lat)
            iPa1 = srcdir2trg[i, 2]           # i + a1
            iMa1 = srcdir2trg[i, 1+L-1]     # i - a1
            iPa2 = srcdir2trg[i, 1+L]       # i + a2
            iMa2 = srcdir2trg[i, 1+L*L-L]  # i - a2

            @simd for k in eachindex(lat)
                kPa1 = srcdir2trg[k, 2]           # k + a1
                kMa1 = srcdir2trg[k, 1+L-1]     # k - a1
                kPa2 = srcdir2trg[k, 1+L]       # k + a2
                kMa2 = srcdir2trg[k, 1+L*L-L]  # k - a2

                temp[1] += 0.25 * weight * measurement.kernel(mc, mc.model,
                               (i, iPa1, iMa1, iPa2, iMa2, k, kPa1, kMa1, kPa2, kMa2), packed_greens, 4)
            end
        end
    end
    return
end


@inline function finalize_temp!(::EachSitePair_B1, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_B1, m)
    push!(m.observable, m.temp[1])
end


###########################
###########################
### magnetic susceptibility Iterators
###########################
###########################

###########################
### lattice iterator for A1 magnetic susceptibility
###########################
"""
    EachSitePair_A1()


"""
struct EachSitePair_A1 <: DirectLatticeIterator end
EachSitePair_A1(::MonteCarloFlavor) = EachSitePair_A1()
output_size(::EachSitePair_A1, l::Lattice) = (1,)

_length(::EachSitePair_A1, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_A1, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachSitePair_A1, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        L = lat.Ls[1]
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i = ix + (iy - 1) * L
            for jx in 1:L, jy in 1:L
                if iseven(ix + jx + iy + jy)
                    j = jx + (jy - 1) * L
                    temp[1] += weight * ((-1)^(ix + jx) + (-1)^(iy + jy)) *
                               measurement.kernel(mc, mc.model, (i, j), packed_greens, 4)
                end
            end
        end

    end
    return
end


@inline function finalize_temp!(::EachSitePair_A1, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_A1, m)
    push!(m.observable, m.temp[1])
end


###########################
###########################
#### Bilinear susceptibility Iterators
###########################
###########################


###########################
### lattice iterator for nematic susceptibility
###########################


"""
    EachDoubleSitePairByDistance()

Creates an iterator template that directly evaluates the nematic susceptibility.

Requires `lattice` to implement `positions` and `lattice_vectors`.
"""
struct EachDoubleSitePairByDistance <: DeferredLatticeIterator end
EachDoubleSitePairByDistance(::MonteCarloFlavor) = EachDoubleSitePairByDistance()

function output_size(::EachDoubleSitePairByDistance, l::Lattice)
    return (1,)
end

function apply!(
    temp::Array, iter::EachDoubleSitePairByDistance, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachDoubleSitePairByDistance, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        srcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L = lat.Ls[1]
        # @inbounds @fastmath for ix in 1:L, iy in 1:L
        #     if isodd(ix+iy)
        #         Δi =(-1)^ix -(-1)^iy
        #         i=ix+(iy-1)*L
        #         for jx in 1:L, jy in 1:L
        #             if isodd(jx+jy)
        #                 Δj=(-1)^jx -(-1)^jy
        #                 j=jx+(jy-1)*L

        #                 tcc=zero(ComplexF64);
        #                 for k in eachindex(lat)  
        #                     kPi = srcdir2trg[k, i]
        #                     @simd for l in eachindex(lat)
        #                         lPj = srcdir2trg[l, j]
        #                         tcc += measurement.kernel(mc, mc.model, (k, l, kPi, lPj), packed_greens, 4)
        #                     end
        #                 end

        #                 temp[1] += weight * Δi*Δj*tcc
        #             end
        #         end
        #     end
        # end
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            if isodd(ix + iy)
                Δi = (-1)^ix - (-1)^iy
                i = ix + (iy - 1) * L
                k = div(N, 2)
                kPi = srcdir2trg[k, i]
                for jx in 1:L, jy in 1:L
                    if isodd(jx + jy)
                        Δj = (-1)^jx - (-1)^jy
                        j = jx + (jy - 1) * L

                        tcc = zero(ComplexF64)
                        @simd for l in eachindex(lat)
                            lPj = srcdir2trg[l, j]
                            tcc += measurement.kernel(mc, mc.model, (k, l, kPi, lPj), packed_greens, 4)
                        end


                        temp[1] += weight * Δi * Δj * tcc * N
                    end
                end
            end
        end
    end
    return
end


@inline function finalize_temp!(::EachDoubleSitePairByDistance, m, mc)
    m.temp ./= (length(lattice(mc))^4)  # *1/N^4
end
@inline function commit!(::EachDoubleSitePairByDistance, m)
    push!(m.observable, m.temp[1])
end



###########################
### lattice iterator for A1_Q1Q2 susceptibility
###########################
"""
EachDoubleSitePairByDistance_Q1Q2()

Creates an iterator template that directly evaluate the A₁` bilinear susceptibility.

Requires `lattice` to implement `positions` and `lattice_vectors`.
"""
struct EachDoubleSitePairByDistance_Q1Q2 <: DeferredLatticeIterator end
EachDoubleSitePairByDistance_Q1Q2(::MonteCarloFlavor) = EachDoubleSitePairByDistance_Q1Q2()

function output_size(::EachDoubleSitePairByDistance_Q1Q2, l::Lattice)
    return (1,)
end



function apply!(
    temp::Array, iter::EachDoubleSitePairByDistance_Q1Q2, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachDoubleSitePairByDistance_Q1Q2, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        srcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L = lat.Ls[1]

        # @inbounds @fastmath for ix in 1:L, iy in 1:L
        #         if iseven(ix+iy)
        #             Δi =(-1)^ix +(-1)^iy
        #             i=ix+(iy-1)*L
        #             for jx in 1:L, jy in 1:L
        #                 if iseven(jx+jy)
        #                     Δj=(-1)^jx +(-1)^jy
        #                     j=jx+(jy-1)*L

        #                     tcc=zero(ComplexF64);
        #                     for k in eachindex(lat)  
        #                         kPi = srcdir2trg[k, i]
        #                         ky, kx=fldmod1(k, L)
        #                         @simd for l in eachindex(lat)
        #                             lPj = srcdir2trg[l, j]
        #                             ly, lx=fldmod1(l, L)
        #                             tcc += (-1)^(lx+ly+kx+ky) * measurement.kernel(mc, mc.model, (k, l, kPi, lPj), packed_greens, 4)
        #                         end
        #                     end

        #                     temp[1] += weight * Δi*Δj*tcc
        #                 end
        #             end
        #         end
        # end 

        k = div(N, 2)
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            if iseven(ix + iy)
                Δi = (-1)^ix + (-1)^iy
                i = ix + (iy - 1) * L
                kPi = srcdir2trg[k, i]
                for jx in 1:L, jy in 1:L
                    if iseven(jx + jy)
                        Δj = (-1)^jx + (-1)^jy
                        j = jx + (jy - 1) * L
                        tcc = zero(ComplexF64)
                        @simd for p in eachindex(lat)
                            kPp = srcdir2trg[k, p]        #l=k+p
                            kPpPj = srcdir2trg[kPp, j]  #l+j=k+p+j
                            py, px = fldmod1(p, L)
                            tcc += (-1)^(px + py) * measurement.kernel(mc, mc.model, (k, kPp, kPi, kPpPj), packed_greens, 4)
                        end
                        temp[1] += weight * Δi * Δj * tcc * N
                    end
                end
            end
        end
    end
    return
end


@inline function finalize_temp!(::EachDoubleSitePairByDistance_Q1Q2, m, mc)
    m.temp ./= (length(lattice(mc))^4)  # *1/N^4
end
@inline function commit!(::EachDoubleSitePairByDistance_Q1Q2, m)
    push!(m.observable, m.temp[1])
end




###########################
### lattice iterator for B1`_Q1Q2 susceptibility
###########################

"""
    EachDoubleSitePairByDistance_B1p_Q1Q2()

Creates an iterator template that directly evaluates the B₁` bilinear susceptibility.

Requires `lattice` to implement `positions` and `lattice_vectors`.
"""
struct EachDoubleSitePairByDistance_B1p_Q1Q2 <: DeferredLatticeIterator end
EachDoubleSitePairByDistance_B1p_Q1Q2(::MonteCarloFlavor) = EachDoubleSitePairByDistance_B1p_Q1Q2()

function output_size(::EachDoubleSitePairByDistance_B1p_Q1Q2, l::Lattice)
    return (1,)
end



function apply!(
    temp::Array, iter::EachDoubleSitePairByDistance_B1p_Q1Q2, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachDoubleSitePairByDistance_B1p_Q1Q2, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L = lat.Ls[1]

        k = div(N, 2)
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            if isodd(ix + iy)
                Δi = (-1)^ix - (-1)^iy
                #note that i and j are treated as directions [so is p], 
                #i.e. in principle, we should subtract (ix,iy) = (ix, iy)- (1, 1). Easy to see, that these subtractions cancel.
                i = ix + (iy - 1) * L
                kPi = Bsrcdir2trg[k, i]
                for jx in 1:L, jy in 1:L
                    if isodd(jx + jy)
                        Δj = (-1)^jx - (-1)^jy
                        j = jx + (jy - 1) * L
                        tcc = zero(ComplexF64)
                        @simd for p in eachindex(lat)
                            kPp = Bsrcdir2trg[k, p]        #l=k+p
                            kPpPj = Bsrcdir2trg[kPp, j]  #l+j=k+p+j
                            py, px = fldmod1(p, L)
                            tcc += (-1)^(px + py) * measurement.kernel(mc, mc.model, (k, kPp, kPi, kPpPj), packed_greens, 4)
                        end
                        temp[1] += weight * Δi * Δj * tcc * N
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
#################
### phase stiffness Iterators
#################
###########################

################################################################################
### Bond pair iterator for phase stiffness 
################################################################################

struct PS_EachBondPairByBravaisDistance <: DeferredLatticeIterator
    bond_idxs::Vector{Int}    #bond indices for b^>, i.e. for the `positive` bonds only.
    dir_idxs::Vector{Int}    #corresp. indices into directions ds.
    bond_weights::Vector{Float64}    #corresp. weights dot(shift_dir, ds[dir_idxs])
    bond_dirs::Vector{Vector{Float64}} #corresp. bond vectors b =ds[dir_idxs]
    bond_cts::Vector{Vector{Float64}}   #corresp. centers of the bond c = ds[dir_idxs]/2
    qs_idxs::Vector{Int}    #linear indices into qs and/or the FFT Λ for the three specific momenta below
    qs_vecs::Vector{Vector{Float64}} # the 3 q-vectors: q_longi, q_trans, q_zero
end

function PS_EachBondPairByBravaisDistance(mc::MonteCarloFlavor; shift_dir=[1, 0])

    normalize!(shift_dir)
    ortho = [-shift_dir[2], shift_dir[1]]

    lat = lattice(mc)
    dirs = directions(lat)
    dir_idxs = filter(hopping_directions(lat)) do i
        dot(dirs[i], shift_dir) > 1e-6
    end
    bond_idxs = hopping_directions_to_bond_idx(dirs[dir_idxs], lat)

    bond_weights = [dot(dirs[i], shift_dir) for i in dir_idxs]

    #q -directions, centered by default, i.e. qs[1,1] elements corresponds to
    # q=(0,0), same as in the FFT.
    qs = cached_reciprocal_discretization(lat)::Matrix{Vector{Float64}}
    #determining q_longitudinal
    idx2dist = map(i -> i => [dot(qs[i], ortho), norm(qs[i])], eachindex(qs))
    filter!(t -> abs(t[2][1]) < 1e-6, idx2dist)
    sort!(idx2dist, by=x -> x[2][2])   #searching for the smallest norm(q)
    idx_longi = idx2dist[2][1]            #picking the second element which corresponds to |q|=2π/L, as the first one is q=(0,0)

    #determining q_transversal
    idx2dist = map(i -> i => [dot(qs[i], shift_dir), norm(qs[i])], eachindex(qs))
    filter!(t -> abs(t[2][1]) < 1e-6, idx2dist)
    sort!(idx2dist, by=x -> x[2][2])
    idx_trans = idx2dist[2][1]

    return PS_EachBondPairByBravaisDistance(bond_idxs, dir_idxs, bond_weights, dirs[dir_idxs], dirs[dir_idxs] * 0.5,
        [idx_longi, idx_trans, 1], [qs[idx_longi], qs[idx_trans], qs[1]])
end


function output_size(iter::PS_EachBondPairByBravaisDistance, l::Lattice)
    if iter.bond_idxs isa Colon
        B = length(l.unitcell.bonds)
    else
        B = length(iter.bond_idxs)
    end
    return (l.Ls..., B, B)
end

function apply!(
    temp::Array, iter::PS_EachBondPairByBravaisDistance,
    measurement, mc::DQMC, packed_greens, weight=1.0
)
    @timeit_debug "apply!(::PS_EachBondPairByBravaisDistance, ::$(typeof(measurement.kernel)))" begin
        l = lattice(mc)
        Lx, Ly = l.Ls
        bs = view(l.unitcell.bonds, iter.bond_idxs)
        modx = get!(l, :modcachex, modcachex)::Vector{Int}
        mody = get!(l, :modcachey, modcachey)::Vector{Int}

        @inbounds @fastmath for σ in measurement.flavor_iterator    #iterates over nothing

            for s1y in 1:Ly, s1x in 1:Lx
                for s2y in 1:Ly, s2x in 1:Lx
                    # output "directions"
                    # 1 because I want onsite at index 1
                    dx = modx[1+s2x-s1x+Lx]
                    dy = mody[1+s2y-s1y+Ly]

                    for (i_b1, b1) in enumerate(bs)
                        s1 = _sub2ind(l, (s1x, s1y, from(b1)))
                        x, y = b1.uc_shift
                        t1 = _sub2ind(l, (modx[s1x+x+Lx], mody[s1y+y+Ly], to(b1)))

                        for (i_b2, b2) in enumerate(bs)
                            s2 = _sub2ind(l, (s2x, s2y, from(b2)))
                            x, y = b2.uc_shift
                            t2 = _sub2ind(l, (modx[s2x+x+Lx], mody[s2y+y+Ly], to(b2)))

                            #(s1, t1, s2, t2)= (r′, r′+ b1, r =r′+d, r+b2)
                            temp[1][dx, dy, i_b1, i_b2] += weight * measurement.kernel(
                                mc, mc.model, (s1, t1, s2, t2), packed_greens, σ
                            )
                        end
                    end
                end
            end
        end
    end
    return
end


@inline function prepare!(::PS_EachBondPairByBravaisDistance, m, mc)
    m.temp[1] .= zero(eltype(m.temp[1]))
    m.temp[2] .= zero(eltype(m.temp[2]))
    m.temp[3] .= zero(eltype(m.temp[3]))
end


@inline function finalize_temp!(::PS_EachBondPairByBravaisDistance, m, mc)
    iter = m.lattice_iterator
    m.temp[1] ./= (length(lattice(mc)))  # *1/N

    for (i, ix_b1) in enumerate(iter.bond_idxs), (j, ix_b2) in enumerate(iter.bond_idxs)
        # Fourier transform on Bravais lattice for on pair of bonds (i, j)
        m.temp[2][:, :] .= m.temp[1][:, :, i, j]
        fft!(m.temp[2])

        # weights from directions and bond centers
        weight = iter.bond_weights[i] * iter.bond_weights[j]
        Δcts = iter.bond_cts[i] - iter.bond_cts[j]
        #we save Λ_L in m.temp[3][1]
        #we save Λ_T in m.temp[3][2]
        #we save Λ[q=(0,0)] in m.temp[3][3]        
        for (nq, ix_q) in enumerate(iter.qs_idxs)
            m.temp[3][nq] += m.temp[2][ix_q] * weight *
                             cis(-dot(Δcts, iter.qs_vecs[nq]))
        end
    end
    #we save the phase stiffness ρₛ in the 4th element
    m.temp[3][4] = 0.25 * (m.temp[3][1] - m.temp[3][2])
end





@inline function commit!(::PS_EachBondPairByBravaisDistance, m)
    #pushing the vector with [Λ_L, Λ_T, Λ[q=(0,0)], ρₛ]
    push!(m.observable, m.temp[3])
end


###########################
### lattice iterator for Kₓ onsite kinetic energy measurement, 
### a.k.a. diamagnetic contribution to phase stiffness
###########################
"""
EachWeightedBond()

By default, we compute Kₓ with `shift_dir=[1, 0]`.
"""
struct EachWeightedBond <: DeferredLatticeIterator
    idxs::Vector{Int}
    weights::Vector{Float64}
end
function EachWeightedBond(mc::MonteCarloFlavor; shift_dir=[1, 0])
    #determining the `positive` bonds, i.e. the bonds with bₓ>0.
    lat = lattice(mc)
    dirs = directions(lat)
    #direction indices idxs --- not the indices into bonds!
    idxs = filter(hopping_directions(lat)) do i
        dot(dirs[i], shift_dir) > 1e-6
    end
    weights = [dot(dirs[i], shift_dir) for i in idxs]
    return EachWeightedBond(idxs, weights)
end


function output_size(::EachWeightedBond, l::Lattice)
    return (1,)
end


function apply!(
    temp::Array, iter::EachWeightedBond, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachWeightedBond, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        dir2srctrg = lat[:dir2srctrg]
        N = length(lattice(mc))
        L = lat.Ls[1]

        @inbounds @fastmath for (ni, i) in enumerate(iter.idxs)  #sum over bonds with bₓ>0
            weight = iter.weights[ni]^2

            for (src, trg) in dir2srctrg[i]     #(src, trg) =(r, r+b)

                temp[1] -= weight * measurement.kernel(mc, mc.model, (src, trg), packed_greens, 4)
            end
        end
    end
    return

end


@inline function finalize_temp!(::EachWeightedBond, m, mc)
    m.temp ./= (length(lattice(mc)))  # *1/N
end
@inline function commit!(::EachWeightedBond, m)
    push!(m.observable, m.temp[1])
end

###########################
###########################
### heat capacity Iterators
###########################
###########################

################################################################################
### Bond pair iterator for heat capacity h₂ 
################################################################################

struct EachDistancedBondPairSummed <: DirectLatticeIterator
    bond_idxs::Vector{Int}
end

function EachDistancedBondPairSummed(mc::MonteCarloFlavor)
    full_numb_bonds = length(lattice(mc).unitcell.bonds)
    bond_idxs = [i for i in 1:full_numb_bonds]
    return EachDistancedBondPairSummed(bond_idxs)
end

function output_size(iter::EachDistancedBondPairSummed, l::Lattice)
    return (1,)
end

function apply!(
    temp::Array, iter::EachDistancedBondPairSummed,
    measurement, mc::DQMC, packed_greens, weight=1.0)
    @timeit_debug "apply!(::EachDistancedBondPairSummed, ::$(typeof(measurement.kernel)))" begin
        l = lattice(mc)
        Lx, Ly = l.Ls
        bs = view(l.unitcell.bonds, iter.bond_idxs)
        modx = get!(l, :modcachex, modcachex)::Vector{Int}
        mody = get!(l, :modcachey, modcachey)::Vector{Int}

        @inbounds @fastmath for σ in measurement.flavor_iterator
            #By default, no flavor iteration.
            for s1y in 1:Ly, s1x in 1:Lx
                for s2y in 1:Ly, s2x in 1:Lx

                    for (i, b1) in enumerate(bs)
                        s1 = _sub2ind(l, (s1x, s1y, from(b1)))
                        x, y = b1.uc_shift
                        t1 = _sub2ind(l, (modx[s1x+x+Lx], mody[s1y+y+Ly], to(b1)))

                        for (j, b2) in enumerate(bs)
                            s2 = _sub2ind(l, (s2x, s2y, from(b2)))
                            x, y = b2.uc_shift
                            t2 = _sub2ind(l, (modx[s2x+x+Lx], mody[s2y+y+Ly], to(b2)))

                            temp[1] += weight * measurement.kernel(mc, mc.model, (s1, t1, s2, t2), packed_greens, σ)
                            #(s1, t1, s2, t2)= (i, i+ b1, j, j+b2)
                        end
                    end
                end
            end
        end
    end
    return
end

@inline function finalize_temp!(::EachDistancedBondPairSummed, m, mc)
    #since we assume that we mostly use the time-averaged measurement, we divide by β
    m.temp ./= (mc.parameters.beta * length(lattice(mc))^2)  # *1/N^2 *1/β
end
@inline function commit!(::EachDistancedBondPairSummed, m)
    push!(m.observable, m.temp[1])
end




################################################################################
### Bond and site iterator for heat capacity h₃
################################################################################

struct EachBondEachSiteSummed <: DirectLatticeIterator
    bond_idxs::Vector{Int}
end

function EachBondEachSiteSummed(mc::MonteCarloFlavor)
    full_numb_bonds = length(lattice(mc).unitcell.bonds)
    bond_idxs = [i for i in 1:full_numb_bonds]
    return EachBondEachSiteSummed(bond_idxs)
end

function output_size(iter::EachBondEachSiteSummed, l::Lattice)
    return (1,)
end

function apply!(
    temp::Array, iter::EachBondEachSiteSummed,
    measurement, mc::DQMC, packed_greens, weight=1.0)
    @timeit_debug "apply!(::EachBondEachSiteSummed, ::$(typeof(measurement.kernel)))" begin
        l = lattice(mc)
        Lx, Ly = l.Ls
        bs = view(l.unitcell.bonds, iter.bond_idxs)
        modx = get!(l, :modcachex, modcachex)::Vector{Int}
        mody = get!(l, :modcachey, modcachey)::Vector{Int}

        @inbounds @fastmath for σ in measurement.flavor_iterator
            #By default, no flavor iteration.
            for s1y in 1:Ly, s1x in 1:Lx
                for (i, b1) in enumerate(bs)
                    s1 = _sub2ind(l, (s1x, s1y, from(b1)))
                    x, y = b1.uc_shift
                    t1 = _sub2ind(l, (modx[s1x+x+Lx], mody[s1y+y+Ly], to(b1)))

                    @simd for j in eachindex(l)
                        temp[1] += weight * measurement.kernel(mc, mc.model, (s1, t1, j), packed_greens, σ)
                        #(s1, t1, j)= (i, i+ b, j)
                    end
                end
            end
        end
    end
    return
end

@inline function finalize_temp!(::EachBondEachSiteSummed, m, mc)
    #since we assume that we mostly use the time-averaged measurement, we divide by β
    U = mc.model.U
    β = mc.parameters.beta
    m.temp .*= (-U / (length(lattice(mc))^2 * β))  # * (-U/N^2) *1/β
end
@inline function commit!(::EachBondEachSiteSummed, m)
    push!(m.observable, m.temp[1])
end



################################################################################
### Site pair iterator for heat capacity h₄
################################################################################

struct EachSiteTwiceSummed <: DirectLatticeIterator
end

function output_size(iter::EachSiteTwiceSummed, l::Lattice)
    return (1,)
end

function apply!(
    temp::Array, iter::EachSiteTwiceSummed,
    measurement, mc::DQMC, packed_greens, weight=1.0)
    @timeit_debug "apply!(::EachSiteTwiceSummed, ::$(typeof(measurement.kernel)))" begin
        l = lattice(mc)
        @inbounds @fastmath for σ in measurement.flavor_iterator
            for i in eachindex(l)
                @simd for j in eachindex(l)
                    temp[1] += weight * measurement.kernel(mc, mc.model, (i, j), packed_greens, σ)
                end
            end
        end
    end
    return
end

@inline function finalize_temp!(::EachSiteTwiceSummed, m, mc)
    U = mc.model.U
    if typeof(m.greens_iterator) <: AbstractUnequalTimeGreensIterator
        β = mc.parameters.beta
    else
        β = 1
    end
    m.temp .*= (U^2 / (length(lattice(mc))^2 * β))  # * U^2/N^2 *1/β
end
@inline function commit!(::EachSiteTwiceSummed, m)
    push!(m.observable, m.temp[1])
end


struct EachSiteSummed <: DirectLatticeIterator
end
function output_size(iter::EachSiteSummed, l::Lattice)
    return (1,)
end

function apply!(
    temp::Array, iter::EachSiteSummed,
    measurement, mc::DQMC, packed_greens, weight=1.0)
    @timeit_debug "apply!(::EachSiteSummed, ::$(typeof(measurement.kernel)))" begin
        l = lattice(mc)
        @inbounds @fastmath for σ in measurement.flavor_iterator
            @simd for i in eachindex(l)
                temp[1] += weight * measurement.kernel(mc, mc.model, (i, i), packed_greens, σ)
            end
        end
    end
    return
end

@inline function finalize_temp!(::EachSiteSummed, m, mc)
    U = mc.model.U
    if typeof(m.greens_iterator) <: AbstractUnequalTimeGreensIterator
        β = mc.parameters.beta
    else
        β = 1
    end
    m.temp .*= (U^2 / (length(lattice(mc))^2 * β))  # * U^2/N^2 *1/β
end
@inline function commit!(::EachSiteSummed, m)
    push!(m.observable, m.temp[1])
end

#############################
### order parameter Iterators
#############################


###########################
### lattice iterator for B₁ bilinear order parameter
###########################

struct EachSitePair_B1_OP <: DirectLatticeIterator end
EachSitePair_B1_OP(::MonteCarloFlavor) = EachSitePair_B1_OP()
output_size(::EachSitePair_B1_OP, l::Lattice) = (1,)
_length(::EachSitePair_B1_OP, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_B1_OP, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachSitePair_B1_OP, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L = lat.Ls[1]
        @inbounds @fastmath for px in 1:L, py in 1:L
            if isodd(px + py)
                Δp = (-1)^(px - 1) - (-1)^(py - 1)    #To convert p to direction we effectively compute px, py =(px, py) -(1,1)
                p = px + (py - 1) * L

                @simd for k in eachindex(lat)
                    kPp = Bsrcdir2trg[k, p]           # k + p                                 
                    temp[1] += weight * Δp * measurement.kernel(mc, mc.model, (k, kPp), packed_greens, 4)
                end
            end
        end
    end
    return
end

@inline function finalize_temp!(::EachSitePair_B1_OP, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_B1_OP, m)
    push!(m.observable, real(m.temp[1]))
end
###########################
### lattice iterator for A₁` bilinear order parameter
###########################

struct EachSitePair_A1p_OP <: DirectLatticeIterator end
EachSitePair_A1p_OP(::MonteCarloFlavor) = EachSitePair_A1p_OP()
output_size(::EachSitePair_A1p_OP, l::Lattice) = (1,)
_length(::EachSitePair_A1p_OP, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_A1p_OP, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachSitePair_A1p_OP, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L = lat.Ls[1]
        @inbounds @fastmath for px in 1:L, py in 1:L
            if iseven(px + py)
                Δp = (-1)^(px - 1) + (-1)^(py - 1)    #To convert p to direction we effectively compute px, py =(px, py) -(1,1)
                p = px + (py - 1) * L

                @simd for k in eachindex(lat)
                    kPp = Bsrcdir2trg[k, p]           # k + p        
                    ky, kx = fldmod1(k, L)
                    temp[1] += weight * Δp * (-1)^(kx + ky) * measurement.kernel(mc, mc.model, (k, kPp), packed_greens, 4)
                end
            end
        end
    end
    return
end

@inline function finalize_temp!(::EachSitePair_A1p_OP, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_A1p_OP, m)
    push!(m.observable, real(m.temp[1]))
end

###########################
### lattice iterator for B₁` bilinear order parameter
###########################

struct EachSitePair_B1p_OP <: DirectLatticeIterator end
EachSitePair_B1p_OP(::MonteCarloFlavor) = EachSitePair_B1p_OP()
output_size(::EachSitePair_B1p_OP, l::Lattice) = (1,)
_length(::EachSitePair_B1p_OP, l::Lattice) = 1

function apply!(
    temp::Array, iter::EachSitePair_B1p_OP, measurement, mc::DQMC,
    packed_greens, weight=1.0
)
    @timeit_debug "apply!(::EachSitePair_B1p_OP, ::$(typeof(measurement.kernel)))" begin
        lat = lattice(mc)
        Bsrcdir2trg = lat[:Bravais_srcdir2trg]::Matrix{Int}

        N = length(Bravais(lat))  #in TwoBandModel N=L^2
        L = lat.Ls[1]
        @inbounds @fastmath for px in 1:L, py in 1:L
            if isodd(px + py)
                Δp = (-1)^(px - 1) - (-1)^(py - 1)    #To convert p to direction we effectively compute px, py =(px, py) -(1,1)
                p = px + (py - 1) * L

                @simd for k in eachindex(lat)
                    kPp = Bsrcdir2trg[k, p]           # k + p        
                    ky, kx = fldmod1(k, L)
                    temp[1] += weight * Δp * (-1)^(kx + ky) * measurement.kernel(mc, mc.model, (k, kPp), packed_greens, 4)
                end
            end
        end
    end
    return
end

@inline function finalize_temp!(::EachSitePair_B1p_OP, m, mc)
    m.temp ./= (length(lattice(mc))^2)  # *1/N^2
end
@inline function commit!(::EachSitePair_B1p_OP, m)
    push!(m.observable, real(m.temp[1]))
end






const _all_lattice_iterator_types = [
    EachSiteAndFlavor, EachSite, EachSitePair, EachSiteByDistance,
    OnSite, EachLocalQuadByDistance, EachLocalQuadBySyncedDistance,
    Sum,
    EachBondPairByBravaisDistance,
    EachSitePair_summed, EachDoubleSitePairByDistance,
    EachDoubleSitePairByDistance_Q1Q2, EachSitePair_B1,
    EachSitePair_B1_OP, EachSitePair_A1p_OP, EachSitePair_B1p_OP,
    EachDoubleSitePairByDistance_B1p_Q1Q2, EachWeightedBond,
    PS_EachBondPairByBravaisDistance, EachDistancedBondPairSummed,
    EachBondEachSiteSummed, EachSiteTwiceSummed
]