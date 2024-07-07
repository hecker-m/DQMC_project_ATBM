"""
    extract_ϕ(mc::DQMC, _config::Int)

Returns the recorded ϕ-configuration field from `mc` at the 
recorded instance `_config`
"""
function extract_ϕ(mc::DQMC, _config::Int)
    δτ = mc.parameters.delta_tau
    my_field = mc.field

    my_field_decompressed = MonteCarlo.decompress(my_field,
        mc.recorder.configs[_config])    #brings the bitstream into array form

    if typeof(my_field) <: AbstractDiscreteMBF
        #convert the values 1,..,4 into the corresponding x-values, 
        # respectively, ϕ~x/sqrt(2U)    cf. field definition
        ϕ_field = my_field.x[my_field_decompressed[:, :, :]]
        ϕ_field ./= sqrt(2δτ)
    else
        ϕ_field = my_field_decompressed
    end
    return ϕ_field
end

"""
get_observable_using_ϕ(_dqmcs::Array, key)

Takes the recorded field configurations (ϕ), 
and computes the quantity corresponding to 'key' at each recorded MC step.
It returns a FullBinner. 
Currently, only implemented for evaluation of the spin susceptibility according to 
Mxϕ_ππ_bar * Mxϕ_ππ_bar
"""
function get_observable_using_ϕ(_dqmcs::Array, value_Q)
    n_workers = length(_dqmcs)
    Nϕ = _dqmcs[1].parameters.Nϕ

    binner_array = make_binner_array(value_Q, Nϕ)
    for worker in 1:n_workers
        number_configs = length(_dqmcs[worker].recorder.configs)

        for config = 1:number_configs
            ϕ_field = extract_ϕ(_dqmcs[worker], config)
            calc_observables!(binner_array, _dqmcs[worker], ϕ_field, value_Q)
        end
    end
    return binner_array
end


"""
For different ordering vectors, there is a different number of observables we can compute.
`make_binner_array()` generates an array with the appropriate number of FullBinners.
"""
make_binner_array(::Val{:Qππ}, Nϕ) = [FullBinner(Float64) for i = 1:4]
make_binner_array(::Val{:Q0πQπ0}, Nϕ) = Nϕ == 1 ? [FullBinner(Float64) for i = 1:(12+2+4+10)] : [FullBinner(Float64) for i = 1:16]
make_binner_array(::Val{:Q0πQπ0_offset}, Nϕ) = Nϕ == 1 ? [FullBinner(Float64) for i = 1:(12+3)] : [FullBinner(Float64) for i = 1:16]


"""
`calc_observables!(binner_array::Array, mc ::DQMC, ϕ_field::Array, ::Val{:Qππ})`

computes all relevant observables corresponding to a magnetic ordering vector Q=(π,π), 
i.e. Neel ordering. The returned array contains `FullBinner` in the following placement: 

#1 is the magnetic order parameter ⟨ |̄ϕ| ⟩ \\
#2 is the magnetic structure factor (1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩ \\
#3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩

"""
function calc_observables!(binner_array::Array, mc::DQMC, ϕ_field::Array, ::Val{:Qππ})
    U = mc.model.U
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Nϕ = mc.parameters.Nϕ
    Nτ = mc.parameters.slices
    δτ = mc.parameters.delta_tau
    β = mc.parameters.beta

    ϕQ3 = calc_ϕQ3(mc, ϕ_field)

    ϕQ3_OP = Vector{Float64}(undef, Nτ)
    ϕQ3_bar = Vector{Float64}(undef, Nϕ)
    ΦA1 = Vector{Float64}(undef, Nτ)
    for ℓ = 1:Nτ
        ϕQ3_OP[ℓ] = norm(ϕQ3[:, ℓ])
        ΦA1[ℓ] = sum(ϕQ3[1:Nϕ, ℓ] .^ 2) / (2U)
    end
    for ζ = 1:Nϕ
        ϕQ3_bar[ζ] = mean(ϕQ3[ζ, :])
    end
    ## #1 is the magnetic order parameter ⟨ |̄ϕ| ⟩
    push!(binner_array[1], mean(ϕQ3_OP))
    ## #2 is the magnetic structure factor S_spin = (1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩
    push!(binner_array[2], mean(ΦA1) - Nϕ / (2U * N * δτ))
    ## #3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩
    push!(binner_array[3], 1 / (2U) * (β * sum(ϕQ3_bar[:] .^ 2) - Nϕ / (N)))
    ## #4 for the Binder cumulant, we also need S_spin^{(2)} = … = ⟨ [ϕ(0)⋅ϕ(0)]² ⟩ = ⟨ ΦA1² ⟩
    push!(binner_array[4], mean(ΦA1 .^ 2) - (2 + Nϕ) / (N * U * δτ) * mean(ΦA1) + Nϕ * (2 + Nϕ) / (4 * δτ^2 * N^2 * U^2))

end


"""
`calc_observables!(binner_array::Array, mc ::DQMC, ϕ_field::Array, ::Val{:Q0πQπ0})`

computes all relevant observables corresponding to a magnetic ordering vector Q=(π,π), 
i.e. Neel ordering. The returned array contains `FullBinner` in the following placement: 

#1 is the magnetic order parameter ⟨ |̄ϕ| ⟩ \\
#2 is the magnetic structure factor (1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩ \\
#3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩

"""
function calc_observables!(binner_array::Array, mc::DQMC, ϕ_field::Array, ::Val{:Q0πQπ0})
    U = mc.model.U
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Nϕ = mc.parameters.Nϕ
    Nτ = mc.parameters.slices
    δτ = mc.parameters.delta_tau
    β = mc.parameters.beta

    ϕQ1Q2 = calc_ϕQ1Q2(mc, ϕ_field)
    ΛA1 = calc_ΛA1(mc, ϕ_field)

    ϕQ1Q2_OP = calc_ϕQ1Q2_OP(mc, ϕQ1Q2)
    ΦA1 = calc_ΦA1(mc, ϕQ1Q2)
    ΦB1 = calc_ΦB1(mc, ϕQ1Q2)
    ΦA1p = calc_ΦA1p(mc, ϕQ1Q2)
    ΦB1p = calc_ΦB1p(mc, ϕQ1Q2)

    ϕQ1Q2_bar = Vector{Float64}(undef, 2Nϕ)
    for ζ = 1:2Nϕ
        ϕQ1Q2_bar[ζ] = mean(ϕQ1Q2[ζ, :])
    end

    ΦPA1, ΦPB1, ΦPA1p = calc_ΦP_bars(mc, ϕ_field)

    S_spin = mean(ΦA1) - Nϕ / (N * U * δτ)
    S2_spin = mean(ΦA1 .^ 2) - 2 * (1 + Nϕ) / (N * U * δτ) * mean(ΦA1) + Nϕ * (1 + Nϕ) / (δτ^2 * N^2 * U^2)
    χ_spin = 1 / (2U) * (β * sum(ϕQ1Q2_bar .^ 2) - 2Nϕ / (N))

    ## #1 is the magnetic order parameter ⟨ |̄ϕ| ⟩ /sqrt(2U)
    push!(binner_array[1], mean(ϕQ1Q2_OP) / sqrt(2U))
    ## #2 is the magnetic structure factor S_spin = ΦA1(0) -Nϕ /(N*U*δτ)   =(1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩
    push!(binner_array[2], S_spin)
    ## #3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩
    push!(binner_array[3], χ_spin)
    ## #4 for the Binder cumulant, we also need S_spin^{(2)} = … = ⟨ [ϕ(0)⋅ϕ(0)]² ⟩ = ⟨ ΦA1² ⟩
    push!(binner_array[4], S2_spin)

    S_bil_B1 = mean(ΦB1 .^ 2) - 2 / (N * U * δτ) * mean(ΦA1) + Nϕ / (δτ^2 * N^2 * U^2)
    S2_bil_B1 = mean(ΦB1 .^ 4) - 12 / (N * δτ * U) * mean(ΦB1 .* ΦB1 .* ΦA1) +
                6 * (Nϕ + 4) / (N^2 * δτ^2 * U^2) * mean(ΦB1 .* ΦB1) + 12 / (N^2 * δτ^2 * U^2) * mean(ΦA1 .* ΦA1) -
                12 * (Nϕ + 2) / (N^3 * δτ^3 * U^3) * mean(ΦA1) + 3Nϕ * (Nϕ + 2) / (N^4 * δτ^4 * U^4)
    χ_bil_B1 = β * mean(ΦB1)^2 - 2 / (N * U) * mean(ΦA1) + Nϕ / (δτ * N^2 * U^2)

    ## #5 is the nematic order parameter ⟨ |ΦB1| ⟩
    push!(binner_array[5], mean(abs.(ΦB1)))
    ## #6 is the nematic structure factor S_{nem}^{B₁} = …  = ⟨ ΦB1²(0) ⟩
    push!(binner_array[6], S_bil_B1)
    ## #7 is the nematic susceptibility β*1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦB1(ℓ) ΦB1(ℓ′)⟩  ± …
    push!(binner_array[7], χ_bil_B1)
    ## #8 for the Binder cumulant, we also need  S_{nem}^{(2),B₁} = …  = ⟨ ΦB1⁴ ⟩
    push!(binner_array[8], S2_bil_B1)

    S_bil_A1p = mean(ΦA1p .^ 2) - 2 / (N * U * δτ) * mean(ΦA1) + Nϕ / (δτ^2 * N^2 * U^2)
    S2_bil_A1p = mean(ΦA1p .^ 4) - 12 / (N * δτ * U) * mean(ΦA1p .* ΦA1p .* ΦA1) +
                 6 * (Nϕ + 4) / (N^2 * δτ^2 * U^2) * mean(ΦA1p .* ΦA1p) + 12 / (N^2 * δτ^2 * U^2) * mean(ΦA1 .* ΦA1) -
                 12 * (Nϕ + 2) / (N^3 * δτ^3 * U^3) * mean(ΦA1) + 3Nϕ * (Nϕ + 2) / (N^4 * δτ^4 * U^4)
    χ_bil_A1p = β * mean(ΦA1p)^2 - 2 / (N * U) * mean(ΦA1) + Nϕ / (δτ * N^2 * U^2)

    ## #9 is the A1′ bilinear order parameter ⟨ |ΦA1′| ⟩
    push!(binner_array[9], mean(abs.(ΦA1p)))
    ## #10 is the A1′ bilinear structure factor ⟨ (ΦA1′)²(0) ⟩
    push!(binner_array[10], S_bil_A1p)
    ## #11 is the A1′ bilinear susceptibility 1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦA1′(ℓ) ΦA1′(ℓ′)⟩  ± …
    push!(binner_array[11], χ_bil_A1p)
    ## #12 for the Binder cumulant, we also need  ⟨ (ΦA1′)⁴ ⟩
    push!(binner_array[12], S2_bil_A1p)


    #h4 binner [13] and h4_OnSite binner [14] 
    h4 = U^2 * (mean(ΛA1 .^ 2) - (Nϕ + 2 / N) / (δτ * U) * mean(ΛA1) + Nϕ * (Nϕ + 2 / N) / (4 * δτ^2 * U^2))
    push!(binner_array[13], h4)
    ## #14 is the interaction energy 
    E_pot = -U * (mean(ΛA1) - 1 / (2 * δτ * U))
    push!(binner_array[14], E_pot)

    push!(binner_array[15], mean(ΦA1p))
    push!(binner_array[16], mean(ΦB1))
    push!(binner_array[17], mean(ϕQ1Q2[1, :]))
    push!(binner_array[18], mean(ϕQ1Q2[2, :]))

    SP_bil_B1 = mean(ΦPB1 .^ 2) - 8 / (N^2 * U * δτ) * mean(ΦPA1) + 4 / (δτ^2 * N^3 * U^2)
    SP2_bil_B1 = mean(ΦPB1 .^ 4) - 48 / (N^2 * δτ * U) * mean(ΦPB1 .* ΦPB1 .* ΦPA1) +
                 24 * (1 + 16 / N) / (N^3 * δτ^2 * U^2) * mean(ΦPB1 .* ΦPB1) + 192 / (N^4 * δτ^2 * U^2) * mean(ΦPA1 .* ΦPA1) -
                 192 * (1 + 8 / N) / (N^5 * δτ^3 * U^3) * mean(ΦPA1) + 48 * (1 + 8 / N) / (N^6 * δτ^4 * U^4)

    χP_bil_B1 = β * mean(ΦPB1)^2 - 8 / (N^2 * U) * mean(ΦPA1) + 4 / (δτ * N^3 * U^2)

    ## #5 is the nematic order parameter ⟨ |ΦB1| ⟩
    push!(binner_array[19], mean(abs.(ΦPB1)))
    ## #6 is the nematic structure factor S_{nem}^{B₁} = …  = ⟨ ΦB1²(0) ⟩
    push!(binner_array[20], SP_bil_B1)
    ## #7 is the nematic susceptibility β*1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦB1(ℓ) ΦB1(ℓ′)⟩  ± …
    push!(binner_array[21], χP_bil_B1)
    ## #8 for the Binder cumulant, we also need  S_{nem}^{(2),B₁} = …  = ⟨ ΦB1⁴ ⟩
    push!(binner_array[22], SP2_bil_B1)

    SP_bil_A1p = mean(ΦPA1p .^ 2) - 8 / (N^2 * U * δτ) * mean(ΦPA1) + 4 / (δτ^2 * N^3 * U^2)
    SP2_bil_A1p = mean(ΦPA1p .^ 4) - 48 / (N^2 * δτ * U) * mean(ΦPA1p .* ΦPA1p .* ΦPA1) +
                  24 * (1 + 16 / N) / (N^3 * δτ^2 * U^2) * mean(ΦPA1p .* ΦPA1p) + 192 / (N^4 * δτ^2 * U^2) * mean(ΦPA1 .* ΦPA1) -
                  192 * (1 + 8 / N) / (N^5 * δτ^3 * U^3) * mean(ΦPA1) + 48 * (1 + 8 / N) / (N^6 * δτ^4 * U^4)

    χP_bil_A1p = β * mean(ΦPA1p)^2 - 8 / (N^2 * U) * mean(ΦPA1) + 4 / (δτ * N^3 * U^2)

    ## #9 is the A1′ bilinear order parameter ⟨ |ΦA1′| ⟩
    push!(binner_array[23], mean(abs.(ΦPA1p)))
    ## #10 is the A1′ bilinear structure factor ⟨ (ΦA1′)²(0) ⟩
    push!(binner_array[24], SP_bil_A1p)
    ## #11 is the A1′ bilinear susceptibility 1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦA1′(ℓ) ΦA1′(ℓ′)⟩  ± …
    push!(binner_array[25], χP_bil_A1p)
    ## #12 for the Binder cumulant, we also need  ⟨ (ΦA1′)⁴ ⟩
    push!(binner_array[26], SP2_bil_A1p)

    push!(binner_array[27], mean(ΦPA1p))
    push!(binner_array[28], mean(ΦPB1))

    if Nϕ > 1
        DA1_ΦB1p = Vector{Float64}(undef, Nτ)
        ΦB1p_OP = zero(Float64)
        for ℓ = 1:Nτ
            DA1_ΦB1p[ℓ] = sum(ΦB1p[:, ℓ] .^ 2)
            ΦB1p_OP += norm(ΦB1p[:, ℓ])
        end
        N_B1p = size(ΦB1p, 1)
        ## #13 is the B1′ bilinear order parameter ⟨ |ΦB1′| ⟩
        push!(binner_array[13], ΦB1p_OP / Nτ)
        ## #14 is the B1′ bilinear structure factor ⟨ [ΦB1′(0)]² ⟩
        push!(binner_array[14], mean(DA1_ΦB1p) - (1 + N_B1p) / (N * U * δτ) * mean(ΦA1) + 2N_B1p / (δτ^2 * N^2 * U^2))
        ## #15 is the B1′ bilinear susceptibility 1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦB1′(ℓ) ΦB1′(ℓ′)⟩  ± …
        push!(binner_array[15], β * sum(mean(ΦB1p, dims=2) .^ 2) - (1 + N_B1p) / (N * U) * mean(ΦA1) + 2N_B1p / (δτ * N^2 * U^2))
        ## #16 for the Binder cumulant, we also need  ⟨ [ΦB1′(0)]⁴ ⟩
        push!(binner_array[16], mean(DA1_ΦB1p .^ 2))
    end
end

function calc_observables!(binner_array::Array, mc::DQMC, ϕ_field::Array, ::Val{:Q0πQπ0_offset})
    U = mc.model.U
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Nϕ = mc.parameters.Nϕ
    Nτ = mc.parameters.slices
    δτ = mc.parameters.delta_tau
    β = mc.parameters.beta

    ϕQ1Q2 = calc_ϕQ1Q2(mc, ϕ_field)
    ΛA1 = calc_ΛA1(mc, ϕ_field)
    Ω4A1 = calc_ΩnA1(mc, ϕ_field, 4)

    ϕQ1Q2_OP = calc_ϕQ1Q2_OP(mc, ϕQ1Q2)
    ΦA1 = calc_ΦA1(mc, ϕQ1Q2)
    ΦB1 = calc_ΦB1(mc, ϕQ1Q2)
    ΦA1p = calc_ΦA1p(mc, ϕQ1Q2)
    ΦB1p = calc_ΦB1p(mc, ϕQ1Q2)

    ϕQ1Q2_bar = calc_ϕQ1Q2_bar(mc, ϕQ1Q2)

    S_spin = mean(ΦA1) - Nϕ / (N * U * δτ)
    S2_spin = mean(ΦA1 .^ 2) - 2 * (1 + Nϕ) / (N * U * δτ) * mean(ΦA1) + Nϕ * (1 + Nϕ) / (δτ^2 * N^2 * U^2)
    χ_spin = 1 / (2U) * (β * sum(ϕQ1Q2_bar .^ 2) - 2Nϕ / (N))

    ## #1 is the magnetic order parameter ⟨ |̄ϕ| ⟩ /sqrt(2U)
    push!(binner_array[1], mean(ϕQ1Q2_OP) / sqrt(2U))
    ## #2 is the magnetic structure factor S_spin = ΦA1(0) -Nϕ /(N*U*δτ)   =(1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩
    push!(binner_array[2], S_spin)
    ## #3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩
    push!(binner_array[3], χ_spin)
    ## #4 for the Binder cumulant, we also need S_spin^{(2)} = … = ⟨ [ϕ(0)⋅ϕ(0)]² ⟩ = ⟨ ΦA1² ⟩
    push!(binner_array[4], S2_spin)

    S_bil_B1 = mean(ΦB1 .^ 2) - 2 / (N * U * δτ) * mean(ΦA1) + Nϕ / (δτ^2 * N^2 * U^2)
    S2_bil_B1 = mean(ΦB1 .^ 4) - 12 / (N * δτ * U) * mean(ΦB1 .* ΦB1 .* ΦA1) +
                6 * (Nϕ + 4) / (N^2 * δτ^2 * U^2) * mean(ΦB1 .* ΦB1) + 12 / (N^2 * δτ^2 * U^2) * mean(ΦA1 .* ΦA1) -
                12 * (Nϕ + 2) / (N^3 * δτ^3 * U^3) * mean(ΦA1) + 3Nϕ * (Nϕ + 2) / (N^4 * δτ^4 * U^4)
    χ_bil_B1 = β * mean(ΦB1)^2 - 2 / (N * U) * mean(ΦA1) + Nϕ / (δτ * N^2 * U^2)

    ## #5 is the nematic order parameter ⟨ |ΦB1| ⟩
    push!(binner_array[5], mean(abs.(ΦB1)))
    ## #6 is the nematic structure factor S_{nem}^{B₁} = …  = ⟨ ΦB1²(0) ⟩
    push!(binner_array[6], S_bil_B1)
    ## #7 is the nematic susceptibility β*1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦB1(ℓ) ΦB1(ℓ′)⟩  ± …
    push!(binner_array[7], χ_bil_B1)
    ## #8 for the Binder cumulant, we also need  S_{nem}^{(2),B₁} = …  = ⟨ ΦB1⁴ ⟩
    push!(binner_array[8], S2_bil_B1)

    S_bil_A1p = mean(ΦA1p .^ 2) - 2 / (N * U * δτ) * mean(ΦA1) + Nϕ / (δτ^2 * N^2 * U^2)
    S2_bil_A1p = mean(ΦA1p .^ 4) - 12 / (N * δτ * U) * mean(ΦA1p .* ΦA1p .* ΦA1) +
                 6 * (Nϕ + 4) / (N^2 * δτ^2 * U^2) * mean(ΦA1p .* ΦA1p) + 12 / (N^2 * δτ^2 * U^2) * mean(ΦA1 .* ΦA1) -
                 12 * (Nϕ + 2) / (N^3 * δτ^3 * U^3) * mean(ΦA1) + 3Nϕ * (Nϕ + 2) / (N^4 * δτ^4 * U^4)
    χ_bil_A1p = β * mean(ΦA1p)^2 - 2 / (N * U) * mean(ΦA1) + Nϕ / (δτ * N^2 * U^2)

    ## #9 is the A1′ bilinear order parameter ⟨ |ΦA1′| ⟩
    push!(binner_array[9], mean(abs.(ΦA1p)))
    ## #10 is the A1′ bilinear structure factor ⟨ (ΦA1′)²(0) ⟩
    push!(binner_array[10], S_bil_A1p)
    ## #11 is the A1′ bilinear susceptibility 1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦA1′(ℓ) ΦA1′(ℓ′)⟩  ± …
    push!(binner_array[11], χ_bil_A1p)
    ## #12 for the Binder cumulant, we also need  ⟨ (ΦA1′)⁴ ⟩
    push!(binner_array[12], S2_bil_A1p)

    #h4 binner [13] and h4_OnSite binner [14] 
    h4 = U^2 * (mean(ΛA1 .^ 2) - (Nϕ + 2 / N) / (δτ * U) * mean(ΛA1) + Nϕ * (Nϕ + 2 / N) / (4 * δτ^2 * U^2))
    h4_OS = U^2 / N * (mean(Ω4A1) - 3 / (δτ * U) * mean(ΛA1) + 3 / (4 * δτ^2 * U^2))
    push!(binner_array[13], h4)
    push!(binner_array[14], h4_OS)
    ## #15 is the interaction energy 
    E_pot = -U * (mean(ΛA1) - 1 / (2 * δτ * U))
    push!(binner_array[15], E_pot)

    if Nϕ > 1
        DA1_ΦB1p = Vector{Float64}(undef, Nτ)
        ΦB1p_OP = zero(Float64)
        for ℓ = 1:Nτ
            DA1_ΦB1p[ℓ] = sum(ΦB1p[:, ℓ] .^ 2)
            ΦB1p_OP += norm(ΦB1p[:, ℓ])
        end
        N_B1p = size(ΦB1p, 1)
        ## #13 is the B1′ bilinear order parameter ⟨ |ΦB1′| ⟩
        push!(binner_array[13], ΦB1p_OP / Nτ)
        ## #14 is the B1′ bilinear structure factor ⟨ [ΦB1′(0)]² ⟩
        push!(binner_array[14], mean(DA1_ΦB1p) - (1 + N_B1p) / (N * U * δτ) * mean(ΦA1) + 2N_B1p / (δτ^2 * N^2 * U^2))
        ## #15 is the B1′ bilinear susceptibility 1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦB1′(ℓ) ΦB1′(ℓ′)⟩  ± …
        push!(binner_array[15], β * sum(mean(ΦB1p, dims=2) .^ 2) - (1 + N_B1p) / (N * U) * mean(ΦA1) + 2N_B1p / (δτ * N^2 * U^2))
        ## #16 for the Binder cumulant, we also need  ⟨ [ΦB1′(0)]⁴ ⟩
        push!(binner_array[16], mean(DA1_ΦB1p .^ 2))
    end
end

"""
    compute_g_d(lat::Lattice, d::Vector)

Computes g(`d`)=(1/(L/2)²) ∑_{p ∈ BZ₄} e^(-i*p*d)

We define the Brillouin zone BZ₄ which is one quarter of the original Brioullin zone, 
as being symmetrically around zero, roughly within px,py  ∈ [-π/2, π/2]. 
By definition, the quantity g(`d`) is real.

"""
function compute_g_d(lat::Lattice)
    L = lat.Ls[1]
    @assert iseven(L) "L=$L is not an even number"
    dirs_Bravais = directions(Bravais(lat))

    dict_p_vals = p_vals_symmetrized(lat)
    p_vals = dict_p_vals[1]
    red_fac = length(dict_p_vals)

    gd_s = Vector{Tuple}()
    for (idx, d) in enumerate(dirs_Bravais)

        gd = zero(ComplexF64)
        for (weight, p) in p_vals
            gd += weight * cis(-dot(p, d))
        end
        @assert isapprox(imag(gd), 0.0, atol=1e-15) "Imaginary part of g(d) is not zero! It is $(imag(gd)) for d=$(d)"

        if !isapprox(real(gd), 0.0, atol=1e-14)
            push!(gd_s, (idx, d, real(gd) / (L^2 / red_fac)))
        end
    end


    return gd_s
end

function p_vals_symmetrized(lat::Lattice{2})
    Q0 = [0, 0]
    Q1 = [π, 0]
    Q2 = [0, π]
    Q3 = [π, π]
    p_vals, dict = p_vals_BZ4(lat; Qs=[Q0, Q1, Q2, Q3])
    return dict
end

function p_vals_BZ4(lat::Lattice{2}; Qs=[[0, 0], [π, 0], [0, π], [π, π]])
    b1, b2 = MonteCarlo.reciprocal_vectors(lat)
    L1, L2 = lat.Ls
    p_vals = [ℓ1 / L1 * b1 + ℓ2 / L2 * b2 for ℓ1 in 0:L1-1, ℓ2 in 0:L2-1][:]
    dict = Dict{Int,Vector{Tuple}}()
    for q in 1:length(Qs)
        dict[q] = Vector{Tuple}()
    end

    for p in p_vals
        temp = [[norm(b1 + b2), (1, 1, 1)]]
        for (idx, Q) in enumerate(Qs), n1 in [0, 1], n2 in [0, 1]
            dist = norm(Q + n1 * b1 + n2 * b2 - p)
            if isapprox(dist, temp[1][1])
                push!(temp, [dist, (idx, n1, n2)])
            elseif dist < temp[1][1]
                temp = [[dist, (idx, n1, n2)]]
            end
        end
        deg = 1 / length(temp)
        for _temp in axes(temp, 1)
            idx, n1, n2 = temp[_temp][2]
            push!(dict[idx], (deg, p - Qs[idx] - n1 * b1 - n2 * b2))
        end
    end
    return p_vals, dict
end






# function generate_p_values_of_small_BZ(lat::Lattice, val::Val{4})
#     p_in_BZ4 = MonteCarlo.cached_reciprocal_discretization(lat)[1:L_half, 1:L_half]
#     offset=0.5*p_in_BZ4[L_half, L_half]
#     for i in eachindex(p_in_BZ4)
#         p_in_BZ4[i] -=offset
#     end
#     return p_in_BZ4
# end
"""
    compute_g_d(lat::Lattice, d::Vector)

Computes g(`d`)=(1/(L/2)²) ∑_{p ∈ BZ₄} e^(-i*p*d)

We define the Brillouin zone BZ₄ which is one quarter of the original Brioullin zone, 
as being symmetrically around zero, roughly within px,py  ∈ [-π/2, π/2]. 
By definition, the quantity g(`d`) is real.

"""
function compute_g_d(lat::Lattice, d::Vector)
    L = lat.Ls[1]
    @assert iseven(L) "L=$L is not an even number"
    L_half = div(L, 2)

    p_in_BZ4 = generate_p_values_of_small_BZ(lat, val)

    gd = zero(ComplexF64)
    for i in eachindex(p_in_BZ4)
        gd += cis(-dot(p_in_BZ4[i], d))
    end

    return real(gd) / L_half^2
end


"""
Calculates 𝜱'(A₁)(ℓ) ,𝜱'(B₁)(ℓ), 𝜱'(A`₁)(ℓ), 
returning the three vectors with ℓ=1,…,N_τ
"""
@inline function calc_ΦP_bars(mc::DQMC, ϕ_field::Array)
    Nϕ = size(ϕ_field, 1)
    U = mc.model.U
    N_slices = mc.parameters.slices
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]
    gd_s = get!(lat, :gd_s, compute_g_d)

    ΦP_bar_A1 = zeros(Float64, N_slices)
    ΦP_bar_B1 = zeros(Float64, N_slices)
    ΦP_bar_A1P = zeros(Float64, N_slices)

    @inbounds @fastmath for (d_index, d_vec, g_d) in gd_s
        for slice in 1:N_slices

            temp_ΦP = zero(Float64)
            temp_ΦP_A1P = zero(Float64)

            for k in axes(ϕ_field, 2)
                ky, kx = fldmod1(k, L)
                kPd = Bsrcdir2trg[k, d_index]           # k + d  
                temp_ΦP += sum(ϕ_field[:, kPd, slice] .* ϕ_field[:, k, slice])
                temp_ΦP_A1P += (-1)^(kx + ky) * sum(ϕ_field[:, kPd, slice] .* ϕ_field[:, k, slice])
            end
            ΦP_bar_A1[slice] += g_d * ((-1)^d_vec[1] + (-1)^d_vec[2]) * temp_ΦP
            ΦP_bar_B1[slice] += g_d * ((-1)^d_vec[1] - (-1)^d_vec[2]) * temp_ΦP
            ΦP_bar_A1P[slice] += g_d * ((-1)^d_vec[1] + (-1)^d_vec[2]) * temp_ΦP_A1P
        end
    end
    return ΦP_bar_A1 ./ (2U * N^2), ΦP_bar_B1 ./ (2U * N^2), ΦP_bar_A1P ./ (2U * N^2)
end



"""
`calc_ϕQ1Q2(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(π,0)] (τ) and ϕ[Q=(0,π)] (τ) returning it as an array ζ=1,..,2*Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ1Q2(mc::DQMC, ϕ_field::Array)
    L = mc.model.l.Ls[1]
    N = length(mc.model.l)
    Nϕ = mc.parameters.Nϕ
    Nτ = mc.parameters.slices
    ϕQ1Q2 = zeros(Float64, 2Nϕ, Nτ)
    for slice in 1:Nτ
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i = ix + (iy - 1) * L
            ϕQ1Q2[1:Nϕ, slice] .+= (-1)^(ix) .* ϕ_field[:, i, slice]
            ϕQ1Q2[Nϕ+1:2Nϕ, slice] .+= (-1)^(iy) .* ϕ_field[:, i, slice]
        end
    end
    return ϕQ1Q2 ./ N
end
"""
`calc_ϕQ1(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(π,0)] (τ) returning it as an array ζ=1,..,Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ1(mc::DQMC, ϕ_field::Array)
    L = mc.model.l.Ls[1]
    N = length(mc.model.l)
    Nϕ = size(ϕ_field, 1)
    Nτ = mc.parameters.slices
    ϕQ1 = zeros(Float64, Nϕ, Nτ)
    for slice in 1:Nτ
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i = ix + (iy - 1) * L
            ϕQ1[:, slice] .+= (-1)^(ix) .* ϕ_field[:, i, slice]
        end
    end
    return ϕQ1 ./ N
end
"""
`calc_ϕQ2(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(0,π)] (τ) returning it as an array ζ=1,..,Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ2(mc::DQMC, ϕ_field::Array)
    L = mc.model.l.Ls[1]
    N = length(mc.model.l)
    Nϕ = size(ϕ_field, 1)
    Nτ = mc.parameters.slices
    ϕQ2 = zeros(Float64, Nϕ, Nτ)
    for slice in 1:Nτ
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i = ix + (iy - 1) * L
            ϕQ2[:, slice] .+= (-1)^(iy) .* ϕ_field[:, i, slice]
        end
    end
    return ϕQ2 ./ N
end
"""
`calc_ϕQ3(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(π,π)] (τ) returning it as an array ζ=1,..,Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ3(mc::DQMC, ϕ_field::Array)
    L = mc.model.l.Ls[1]
    N = length(mc.model.l)
    Nϕ = size(ϕ_field, 1)
    Nτ = mc.parameters.slices
    ϕQ3 = zeros(Float64, Nϕ, Nτ)
    for slice in 1:Nτ
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i = ix + (iy - 1) * L
            ϕQ3[:, slice] .+= (-1)^(ix + iy) .* ϕ_field[:, i, slice]
        end
    end
    return ϕQ3 ./ N
end

@inline function calc_ΛA1(mc::DQMC, ϕ_field::Array)
    lat = mc.model.l
    N = length(lat)
    U = mc.model.U
    Nτ = mc.parameters.slices
    ΛA1 = zeros(Float64, Nτ)
    for slice in 1:Nτ
        @inbounds @fastmath for i in eachindex(lat)
            ΛA1[slice] += sum(ϕ_field[:, i, slice] .^ 2)
        end
    end
    return ΛA1 ./ (2U * N)
end

@inline function calc_ΩnA1(mc::DQMC, ϕ_field::Array, n::Int)
    lat = mc.model.l
    N = length(lat)
    U = mc.model.U
    Nτ = mc.parameters.slices
    ΩA1 = zeros(Float64, Nτ)
    for slice in 1:Nτ
        @inbounds @fastmath for i in eachindex(lat)
            ΩA1[slice] += sum(ϕ_field[:, i, slice] .^ n)
        end
    end
    return ΩA1 ./ ((2U)^(n / 2) * N)
end

####################
### bilinears of ϕ=(ϕ_Q₁, ϕ_Q₂)
####################
@inline function calc_ϕQ1Q2_OP(mc::DQMC, ϕQ1Q2::Array)
    Nτ = mc.parameters.slices
    ϕQ1Q2_OP = Vector{Float64}(undef, Nτ)
    for ℓ = 1:Nτ
        ϕQ1Q2_OP[ℓ] = norm(ϕQ1Q2[:, ℓ])
    end
    return ϕQ1Q2_OP
end
@inline function calc_ϕQ1Q2_bar(mc::DQMC, ϕQ1Q2::Array)
    Nϕ = mc.parameters.Nϕ
    ϕQ1Q2_bar = Vector{Float64}(undef, 2Nϕ)
    for ζ = 1:2Nϕ
        ϕQ1Q2_bar[ζ] = mean(ϕQ1Q2[ζ, :])
    end
    return ϕQ1Q2_bar
end

@inline function calc_ΦA1(mc::DQMC, ϕQ1Q2::Array)
    Nϕ = mc.parameters.Nϕ
    U = mc.model.U
    Nτ = mc.parameters.slices
    ΦA1 = Vector{Float64}(undef, Nτ)
    for ℓ = 1:Nτ
        ΦA1[ℓ] = sum(ϕQ1Q2[1:Nϕ, ℓ] .^ 2 + ϕQ1Q2[Nϕ+1:2Nϕ, ℓ] .^ 2) / (2U)
    end
    return ΦA1
end
@inline function calc_ΦB1(mc::DQMC, ϕQ1Q2::Array)
    Nϕ = mc.parameters.Nϕ
    U = mc.model.U
    Nτ = mc.parameters.slices
    ΦB1 = Vector{Float64}(undef, Nτ)
    for ℓ = 1:Nτ
        ΦB1[ℓ] = sum(ϕQ1Q2[1:Nϕ, ℓ] .^ 2 - ϕQ1Q2[Nϕ+1:2Nϕ, ℓ] .^ 2) / (2U)
    end
    return ΦB1
end
@inline function calc_ΦA1p(mc::DQMC, ϕQ1Q2::Array)
    Nϕ = mc.parameters.Nϕ
    U = mc.model.U
    Nτ = mc.parameters.slices
    ΦA1p = Vector{Float64}(undef, Nτ)
    for ℓ = 1:Nτ
        ΦA1p[ℓ] = 2 * sum(ϕQ1Q2[1:Nϕ, ℓ] .* ϕQ1Q2[Nϕ+1:2Nϕ, ℓ]) / (2U)
    end
    return ΦA1p
end
@inline function calc_ΦB1p(mc::DQMC, ϕQ1Q2::Array)
    Nϕ = mc.parameters.Nϕ
    U = mc.model.U
    Nτ = mc.parameters.slices
    if Nϕ == 1
        # println("Φ_B1` cannot be computed from ϕ for Nϕ=$(Nϕ)")
        return nothing
    end
    if Nϕ == 2
        ΦB1p = Array{Float64}(undef, 1, Nτ)
    elseif Nϕ == 3
        ΦB1p = Array{Float64}(undef, 3, Nτ)
    end

    for ℓ = 1:Nτ
        if Nϕ == 2
            ΦB1p[1, ℓ] = 2 * (ϕQ1Q2[1, ℓ] * ϕQ1Q2[Nϕ+2, ℓ] - ϕQ1Q2[2, ℓ] * ϕQ1Q2[Nϕ+1, ℓ]) / (2U)
        elseif Nϕ == 3
            ΦB1p[1, ℓ] = 2 * (ϕQ1Q2[2, ℓ] * ϕQ1Q2[Nϕ+3, ℓ] - ϕQ1Q2[3, ℓ] * ϕQ1Q2[Nϕ+2, ℓ]) / (2U)
            ΦB1p[2, ℓ] = 2 * (ϕQ1Q2[3, ℓ] * ϕQ1Q2[Nϕ+1, ℓ] - ϕQ1Q2[1, ℓ] * ϕQ1Q2[Nϕ+3, ℓ]) / (2U)
            ΦB1p[3, ℓ] = 2 * (ϕQ1Q2[1, ℓ] * ϕQ1Q2[Nϕ+2, ℓ] - ϕQ1Q2[2, ℓ] * ϕQ1Q2[Nϕ+1, ℓ]) / (2U)
        end
    end
    return ΦB1p
end

######################################################
#### older functions
######################################################
"""
Calculates ̄M^(xϕ)[Q=(π,π)], directly returning the scalar value
"""
function calc_Mxϕ_ππ_bar(mc::DQMC, ϕ_field::Array)
    U = mc.model.U
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Mzϕ_ππ_bar = zero(Float64)

    @inbounds @fastmath for ix in 1:L, iy in 1:L
        i = ix + (iy - 1) * L
        Mzϕ_ππ_bar += (-1)^(ix + iy) * mean(ϕ_field[1, i, :])
    end
    Mzϕ_ππ_bar /= (sqrt(2U) * N)
    return Mzϕ_ππ_bar
end


"""
Calculates Φ^(B₁)_ϕ (ℓ), returning the vector ℓ=1,…,N_τ
"""
function calc_ΦB1_bar(mc::DQMC, ϕ_field::Array)
    Nϕ = size(ϕ_field, 1)
    U = mc.model.U
    N_slices = mc.parameters.slices
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]

    ϕB1_bar = zeros(Float64, N_slices)
    # for slice in 1:N_slices
    #     for i in 1:N
    #         iy, ix=fldmod1(i, L)                                    
    #         if isodd(ix+iy)
    #             Δmi=((-1)^(ix-1)-(-1)^(iy-1))
    #             temp0=zero(Float64)
    #             for k in 1:N
    #                 kPi=Bsrcdir2trg[k, i]           # k + i  
    #                 temp0 += sum(ϕ_field[: ,kPi ,slice] .*ϕ_field[:, k,slice])
    #             end
    #             ϕB1_bar[slice] +=Δmi * temp0
    #         end
    #     end
    # end
    for slice in 1:N_slices
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            if isodd(ix + iy)
                i = ix + (iy - 1) * L
                Δmi = (-1)^(ix - 1) - (-1)^(iy - 1)
                temp0 = zero(Float64)
                @simd for k in eachindex(lat)
                    kPi = Bsrcdir2trg[k, i]           # k + i  
                    temp0 += sum(ϕ_field[:, kPi, slice] .* ϕ_field[:, k, slice])
                end
                ϕB1_bar[slice] += Δmi * temp0
            end
        end
    end
    ϕB1_bar ./= (2U * N^2)
    return ϕB1_bar
end


"""
Calculates ̄Φ⁰_ϕ (ℓ), returning the vector ℓ=1,…,N_τ
"""
function calc_Φ0_bar(mc::DQMC, ϕ_field::Array)
    Nϕ = size(ϕ_field, 1)
    U = mc.model.U
    N_slices = mc.parameters.slices
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]

    ϕ0_bar = zeros(Float64, N_slices)
    for slice in 1:N_slices
        for i in 1:N
            iy, ix = fldmod1(i, L)
            if iseven(ix + iy)
                Δpi = ((-1)^(ix - 1) + (-1)^(iy - 1))
                temp0 = zero(Float64)
                for k in 1:N
                    kPi = Bsrcdir2trg[k, i]           # k + i  
                    temp0 += sum(ϕ_field[:, kPi, slice] .* ϕ_field[:, k, slice])
                end
                ϕ0_bar[slice] += Δpi * temp0
            end
        end
    end
    return ϕ0_bar ./ (2U * N^2)
end

"""
Calculates ̄Φ^(A₁')_ϕ (ℓ), returning the vector ℓ=1,…,N_τ
"""
function calc_ΦA1p_bar(mc::DQMC, ϕ_field::Array)
    Nϕ = size(ϕ_field, 1)
    U = mc.model.U
    N_slices = mc.parameters.slices
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]

    ϕA1p_bar = zeros(Float64, N_slices)
    for slice in 1:N_slices
        for i in 1:N
            iy, ix = fldmod1(i, L)
            if iseven(ix + iy)
                Δpi = ((-1)^(ix - 1) + (-1)^(iy - 1))
                temp0 = zero(Float64)
                for k in 1:N
                    ky, kx = fldmod1(k, L)
                    kPi = Bsrcdir2trg[k, i]           # k + i  
                    temp0 += (-1)^(kx + ky) * sum(ϕ_field[:, kPi, slice] .* ϕ_field[:, k, slice])
                end
                ϕA1p_bar[slice] += Δpi * temp0
            end
        end
    end
    return ϕA1p_bar ./ (2U * N^2)
end

"""
Calculates ̄Φ^(B₁')_ϕ (ℓ), returning the matrix, η=1,…,Nϕ, ℓ=1,…,N_τ
"""
function calc_ΦB1p_bar(mc::DQMC, ϕ_field::Array)
    Nϕ = size(ϕ_field, 1)
    U = mc.model.U
    N_slices = mc.parameters.slices
    lat = mc.model.l
    L = lat.Ls[1]
    N = length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]
    if Nϕ == 1
        println("Φ^(B₁') is not defined for the Ising case with Nϕ=1")
        return nothing
    end
    if Nϕ == 2
        ϕB1p_bar = zeros(Float64, 1, N_slices)
        for slice in 1:N_slices
            for i in 1:N
                iy, ix = fldmod1(i, L)
                if isodd(ix + iy)
                    Δmi = ((-1)^(ix - 1) - (-1)^(iy - 1))
                    temp0 = zero(Float64)
                    for k in 1:N
                        ky, kx = fldmod1(k, L)
                        kPi = Bsrcdir2trg[k, i]           # k + i  
                        temp0 += (-1)^(kx + ky) * (ϕ_field[1, kPi, slice] * ϕ_field[2, k, slice] - ϕ_field[2, kPi, slice] * ϕ_field[1, k, slice])
                    end
                    ϕB1p_bar[1, slice] += Δmi * temp0
                end
            end
        end
        return ϕB1p_bar ./ (2U * N^2)
    else
        ϕB1p_bar = zeros(Float64, 3, N_slices)
        for slice in 1:N_slices
            for i in 1:N
                iy, ix = fldmod1(i, L)
                if isodd(ix + iy)
                    Δmi = ((-1)^(ix - 1) - (-1)^(iy - 1))
                    temp0 = zeros(Float64, 3)
                    for k in 1:N
                        ky, kx = fldmod1(k, L)
                        kPi = Bsrcdir2trg[k, i]           # k + i  
                        temp0[1] += (-1)^(kx + ky) * (ϕ_field[2, kPi, slice] * ϕ_field[3, k, slice] - ϕ_field[3, kPi, slice] * ϕ_field[2, k, slice])
                        temp0[2] += (-1)^(kx + ky) * (ϕ_field[3, kPi, slice] * ϕ_field[1, k, slice] - ϕ_field[1, kPi, slice] * ϕ_field[3, k, slice])
                        temp0[3] += (-1)^(kx + ky) * (ϕ_field[1, kPi, slice] * ϕ_field[2, k, slice] - ϕ_field[2, kPi, slice] * ϕ_field[1, k, slice])
                    end
                    ϕB1p_bar[:, slice] .+= Δmi .* temp0[:]
                end
            end
        end
        return ϕB1p_bar ./ (2U * N^2)
    end
end