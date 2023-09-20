"""
get_observable_using_ϕ(_dqmcs::Array, key)

Takes the recorded field configurations (ϕ), 
and computes the quantity corresponding to 'key' at each recorded MC step.
It returns a FullBinner. 
Currently, only implemented for evaluation of the spin susceptibility according to 
Mxϕ_ππ_bar * Mxϕ_ππ_bar
"""
function get_observable_using_ϕ(_dqmcs::Array, value_Q)
    n_workers=length(_dqmcs)
    values=FullBinner(Float64)
    Nϕ=_dqmcs[1].parameters.Nϕ
    N=length(_dqmcs[1].model.l)
    U=_dqmcs[1].model.U
    Nτ=_dqmcs[1].parameters.slices
    δτ=_dqmcs[1].parameters.delta_tau
    β=_dqmcs[1].parameters.beta

    binner_array=make_binner_array(value_Q, Nϕ)
    for worker in 1:n_workers
        my_field=_dqmcs[worker].field
        number_configs=length(_dqmcs[worker].recorder.configs)

        for config=1:number_configs
            my_field_decompressed=MonteCarlo.decompress(my_field, 
                _dqmcs[worker].recorder.configs[config])    #brings the bitstream into array form
            ϕ_field=my_field.x[my_field_decompressed[:,:,:]]  
            ϕ_field ./=sqrt(2δτ)  
            #convert the values 1,..,4 into the corresponding x-values, 
            # respectively, ϕ~x/sqrt(2U)    cf. field definition

            calc_observables!(binner_array, _dqmcs[worker], ϕ_field, value_Q)
          
        end
    end
    return binner_array
end


"""
For different ordering vectors, there is a different number of observables we can compute.
`make_binner_array()` generates an array with the appropriate number of FullBinners.
"""
make_binner_array(::Val{:Qππ}, Nϕ)=[FullBinner(Float64) for i=1:4]
make_binner_array(::Val{:Q0πQπ0}, Nϕ)= Nϕ==1 ? [FullBinner(Float64) for i=1:12] : [FullBinner(Float64) for i=1:16]


"""
`calc_observables!(binner_array::Array, mc ::DQMC, ϕ_field::Array, ::Val{:Qππ})`

computes all relevant observables corresponding to a magnetic ordering vector Q=(π,π), 
i.e. Neel ordering. The returned array contains `FullBinner` in the following placement: 

#1 is the magnetic order parameter ⟨ |̄ϕ| ⟩ \\
#2 is the magnetic structure factor (1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩ \\
#3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩

"""
function calc_observables!(binner_array::Array, mc ::DQMC, ϕ_field::Array, ::Val{:Qππ})
    U=mc.model.U
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Nϕ=mc.parameters.Nϕ
    Nτ=mc.parameters.slices
    δτ=mc.parameters.delta_tau
    β=mc.parameters.beta

    ϕQ3= calc_ϕQ3(mc, ϕ_field)

    ϕQ3_OP=Vector{Float64}(undef, Nτ)
    ϕQ3_bar=Vector{Float64}(undef, Nϕ)
    ΦA1=Vector{Float64}(undef, Nτ)
    for ℓ=1:Nτ
        ϕQ3_OP[ℓ]=norm(ϕQ3[:,ℓ])
        ΦA1[ℓ]=sum(ϕQ3[1:Nϕ,ℓ] .^2 )/(2U)
    end
    for ζ=1:Nϕ
        ϕQ3_bar[ζ]=mean(ϕQ3[ζ,:])
    end
    ## #1 is the magnetic order parameter ⟨ |̄ϕ| ⟩
    push!(binner_array[1], mean(ϕQ3_OP))
    ## #2 is the magnetic structure factor S_spin = (1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩
    push!(binner_array[2], mean(ΦA1) -Nϕ/(2U*N*δτ))
    ## #3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩
    push!(binner_array[3], 1/(2U) *(β*sum(ϕQ3_bar[:] .^2) - Nϕ /(N)))
    ## #4 for the Binder cumulant, we also need S_spin^{(2)} = … = ⟨ [ϕ(0)⋅ϕ(0)]² ⟩ = ⟨ ΦA1² ⟩
    push!(binner_array[4], mean(ΦA1 .^2)  - (2+Nϕ)/(N*U*δτ) * mean(ΦA1) + Nϕ*(2+Nϕ)/(4*δτ^2*N^2 *U^2)  ) 
  
end


"""
`calc_observables!(binner_array::Array, mc ::DQMC, ϕ_field::Array, ::Val{:Q0πQπ0})`

computes all relevant observables corresponding to a magnetic ordering vector Q=(π,π), 
i.e. Neel ordering. The returned array contains `FullBinner` in the following placement: 

#1 is the magnetic order parameter ⟨ |̄ϕ| ⟩ \\
#2 is the magnetic structure factor (1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩ \\
#3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩

"""
function calc_observables!(binner_array::Array, mc ::DQMC, ϕ_field::Array, ::Val{:Q0πQπ0})
    U=mc.model.U
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Nϕ=mc.parameters.Nϕ
    Nτ=mc.parameters.slices
    δτ=mc.parameters.delta_tau
    β=mc.parameters.beta

    ϕQ1Q2= calc_ϕQ1Q2(mc, ϕ_field)

    ϕQ1Q2_OP=Vector{Float64}(undef, Nτ)
    ΦA1=Vector{Float64}(undef, Nτ)
    ΦB1=Vector{Float64}(undef, Nτ)
    ΦA1p=Vector{Float64}(undef, Nτ)
    if Nϕ==2
        ΦB1p=Array{Float64}(undef,1, Nτ)
    elseif Nϕ==3
        ΦB1p=Array{Float64}(undef,3, Nτ)
    end
    for ℓ=1:Nτ
        ϕQ1Q2_OP[ℓ]=norm(ϕQ1Q2[:,ℓ])
        ΦA1[ℓ]=sum(ϕQ1Q2[1:Nϕ,ℓ] .^2  +  ϕQ1Q2[Nϕ+1:2Nϕ,ℓ] .^2)/(2U)
        ΦB1[ℓ]=sum(ϕQ1Q2[1:Nϕ,ℓ] .^2  -  ϕQ1Q2[Nϕ+1:2Nϕ,ℓ] .^2)/(2U)
        ΦA1p[ℓ]=2*sum(ϕQ1Q2[1:Nϕ,ℓ] .*  ϕQ1Q2[Nϕ+1:2Nϕ,ℓ])/(2U)
        if Nϕ==2
            ΦB1p[1,ℓ]=2*(ϕQ1Q2[1,ℓ]*ϕQ1Q2[Nϕ+2,ℓ]-ϕQ1Q2[2,ℓ]*ϕQ1Q2[Nϕ+1,ℓ])/(2U)
        elseif Nϕ==3
            ΦB1p[1,ℓ]=2*(ϕQ1Q2[2,ℓ]*ϕQ1Q2[Nϕ+3,ℓ]-ϕQ1Q2[3,ℓ]*ϕQ1Q2[Nϕ+2,ℓ])/(2U)
            ΦB1p[2,ℓ]=2*(ϕQ1Q2[3,ℓ]*ϕQ1Q2[Nϕ+1,ℓ]-ϕQ1Q2[1,ℓ]*ϕQ1Q2[Nϕ+3,ℓ])/(2U)
            ΦB1p[3,ℓ]=2*(ϕQ1Q2[1,ℓ]*ϕQ1Q2[Nϕ+2,ℓ]-ϕQ1Q2[2,ℓ]*ϕQ1Q2[Nϕ+1,ℓ])/(2U)
        end
    end
    ϕQ1Q2_bar=Vector{Float64}(undef, 2Nϕ)
    for ζ=1:2Nϕ
        ϕQ1Q2_bar[ζ]=mean(ϕQ1Q2[ζ,:])
    end
    ## #1 is the magnetic order parameter ⟨ |̄ϕ| ⟩ /sqrt(2U)
    push!(binner_array[1], mean(ϕQ1Q2_OP)/sqrt(2U))
    ## #2 is the magnetic structure factor S_spin = ΦA1(0) -Nϕ /(N*U*δτ)   =(1/(2U)) *∑_ζ ⟨(1/Nτ)∑_ℓ ϕ²(ℓ) -1/(N*δτ)⟩
    push!(binner_array[2], mean(ΦA1) -Nϕ/(N*U*δτ) )
    ## #3 is the magnetic susceptibility (1/(2U)) *∑_ζ ⟨β  ̄ϕ²  -1/(N)⟩
    push!(binner_array[3], 1/(2U) *(β* sum(ϕQ1Q2_bar .^2) - 2Nϕ /(N)))
    ## #4 for the Binder cumulant, we also need S_spin^{(2)} = … = ⟨ [ϕ(0)⋅ϕ(0)]² ⟩ = ⟨ ΦA1² ⟩
    push!(binner_array[4], mean(ΦA1 .^2)  - 2*(1+Nϕ)/(N*U*δτ) * mean(ΦA1) + Nϕ*(1+Nϕ)/(δτ^2*N^2 *U^2))

    ## #5 is the nematic order parameter ⟨ |ΦB1| ⟩
    push!(binner_array[5], mean(abs.(ΦB1)))    
    ## #6 is the nematic structure factor S_{nem}^{B₁} = …  = ⟨ ΦB1²(0) ⟩
    push!(binner_array[6], mean(ΦB1 .^2) - 2/(N*U*δτ) * mean(ΦA1) + Nϕ/(δτ^2*N^2 *U^2))
    ## #7 is the nematic susceptibility β*1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦB1(ℓ) ΦB1(ℓ′)⟩  ± …
    push!(binner_array[7], β*mean(ΦB1)^2 - 2/(N*U) * mean(ΦA1) + Nϕ/(δτ*N^2 *U^2))    
    ## #8 for the Binder cumulant, we also need  ⟨ ΦB1⁴ ⟩
    push!(binner_array[8], mean(ΦB1 .^4)  )

    ## #9 is the A1′ bilinear order parameter ⟨ |ΦA1′| ⟩
    push!(binner_array[9], mean(abs.(ΦA1p)))    
    ## #10 is the A1′ bilinear structure factor ⟨ (ΦA1′)²(0) ⟩
    push!(binner_array[10], mean(ΦA1p .^2) - 2/(N*U*δτ) * mean(ΦA1) + Nϕ/(δτ^2*N^2 *U^2))
    ## #11 is the A1′ bilinear susceptibility 1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦA1′(ℓ) ΦA1′(ℓ′)⟩  ± …
    push!(binner_array[11], β*mean(ΦA1p)^2 - 2/(N*U) * mean(ΦA1) + Nϕ/(δτ*N^2 *U^2)) 
    ## #12 for the Binder cumulant, we also need  ⟨ (ΦA1′)⁴ ⟩
    push!(binner_array[12], mean(ΦA1p .^4)  )
    if Nϕ>1
        DA1_ΦB1p=Vector{Float64}(undef, Nτ)
        ΦB1p_OP=zero(Float64)
        for ℓ=1:Nτ
            DA1_ΦB1p[ℓ]=sum(ΦB1p[:,ℓ].^2)
            ΦB1p_OP+=norm(ΦB1p[:,ℓ])
        end
        N_B1p=size(ΦB1p, 1)
        ## #13 is the B1′ bilinear order parameter ⟨ |ΦB1′| ⟩
        push!(binner_array[13], ΦB1p_OP/Nτ)    
        ## #14 is the B1′ bilinear structure factor ⟨ [ΦB1′(0)]² ⟩
        push!(binner_array[14], mean(DA1_ΦB1p) - (1+N_B1p)/(N*U*δτ) * mean(ΦA1) + 2N_B1p/(δτ^2*N^2 *U^2) )
        ## #15 is the B1′ bilinear susceptibility 1/(Nτ²)∑_{ℓ,ℓ′} ⟨ ΦB1′(ℓ) ΦB1′(ℓ′)⟩  ± …
        push!(binner_array[15], β*sum(mean(ΦB1p, dims=2) .^2) - (1+N_B1p)/(N*U) * mean(ΦA1) +2N_B1p/(δτ*N^2 *U^2)) 
        ## #16 for the Binder cumulant, we also need  ⟨ [ΦB1′(0)]⁴ ⟩
        push!(binner_array[16], mean(DA1_ΦB1p .^2)  )
    end
end



"""
`calc_ϕQ1Q2(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(π,0)] (τ) and ϕ[Q=(0,π)] (τ) returning it as an array ζ=1,..,2*Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ1Q2(mc ::DQMC, ϕ_field::Array)
    L=mc.model.l.Ls[1]
    N=length(mc.model.l)
    Nϕ=mc.parameters.Nϕ
    Nτ=mc.parameters.slices
    ϕQ1Q2= zeros(Float64, 2Nϕ, Nτ)
    for slice in 1:Nτ    
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i=ix+(iy-1)*L
            ϕQ1Q2[1:Nϕ, slice] .+= (-1)^(ix) .* ϕ_field[: ,i ,slice]
            ϕQ1Q2[Nϕ+1:2Nϕ, slice] .+= (-1)^(iy) .* ϕ_field[: ,i ,slice]
        end
    end
    return ϕQ1Q2 ./N
end
"""
`calc_ϕQ1(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(π,0)] (τ) returning it as an array ζ=1,..,Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ1(mc ::DQMC, ϕ_field::Array)
    L=mc.model.l.Ls[1]
    N=length(mc.model.l)
    Nϕ=size(ϕ_field, 1)
    Nτ=mc.parameters.slices
    ϕQ1= zeros(Float64, Nϕ, Nτ)
    for slice in 1:Nτ    
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i=ix+(iy-1)*L
            ϕQ1[:, slice] .+= (-1)^(ix) .* ϕ_field[: ,i ,slice]
        end
    end
    return ϕQ1 ./N
end
"""
`calc_ϕQ2(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(0,π)] (τ) returning it as an array ζ=1,..,Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ2(mc ::DQMC, ϕ_field::Array)
    L=mc.model.l.Ls[1]
    N=length(mc.model.l)
    Nϕ=size(ϕ_field, 1)
    Nτ=mc.parameters.slices
    ϕQ2= zeros(Float64, Nϕ, Nτ)
    for slice in 1:Nτ    
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i=ix+(iy-1)*L
            ϕQ2[:, slice] .+= (-1)^(iy) .* ϕ_field[: ,i ,slice]
        end
    end
    return ϕQ2 ./N
end
"""
`calc_ϕQ3(mc ::DQMC, ϕ_field::Array)` computes the Fourier transform 
ϕ[Q=(π,π)] (τ) returning it as an array ζ=1,..,Nϕ and ℓ=1,..,Nτ 
"""
@inline function calc_ϕQ3(mc ::DQMC, ϕ_field::Array)
    L=mc.model.l.Ls[1]
    N=length(mc.model.l)
    Nϕ=size(ϕ_field, 1)
    Nτ=mc.parameters.slices
    ϕQ3= zeros(Float64, Nϕ, Nτ)
    for slice in 1:Nτ    
        @inbounds @fastmath for ix in 1:L, iy in 1:L
            i=ix+(iy-1)*L
            ϕQ3[:, slice] .+= (-1)^(ix+iy) .* ϕ_field[: ,i ,slice]
        end
    end
    return ϕQ3 ./N
end
"""
Calculates ̄M^(xϕ)[Q=(π,π)], directly returning the scalar value
"""
function calc_Mxϕ_ππ_bar(mc ::DQMC, ϕ_field::Array)
    U=mc.model.U
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Mzϕ_ππ_bar= zero(Float64)

    @inbounds @fastmath for ix in 1:L, iy in 1:L
        i=ix+(iy-1)*L
        Mzϕ_ππ_bar += (-1)^(ix+iy) * mean(ϕ_field[1 ,i ,:])
    end
    Mzϕ_ππ_bar /= (sqrt(2U)*N)
    return Mzϕ_ππ_bar 
end


"""
Calculates Φ^(B₁)_ϕ (ℓ), returning the vector ℓ=1,…,N_τ
"""
function calc_ΦB1_bar(mc ::DQMC, ϕ_field::Array)
    Nϕ=size(ϕ_field, 1)
    U=mc.model.U
    N_slices=mc.parameters.slices
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]

    ϕB1_bar= zeros(Float64, N_slices)
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
            if isodd(ix+iy)
                i=ix+(iy-1)*L
                Δmi=(-1)^(ix-1)-(-1)^(iy-1)
                temp0=zero(Float64)
                @simd for k in eachindex(lat)
                    kPi=Bsrcdir2trg[k, i]           # k + i  
                    temp0 += sum(ϕ_field[: ,kPi ,slice] .*ϕ_field[:, k,slice])
                end
                ϕB1_bar[slice] += Δmi * temp0
            end
        end
    end
    ϕB1_bar ./= (2U*N^2)
    return ϕB1_bar 
end


"""
Calculates ̄Φ⁰_ϕ (ℓ), returning the vector ℓ=1,…,N_τ
"""
function calc_Φ0_bar(mc ::DQMC, ϕ_field::Array)
    Nϕ=size(ϕ_field, 1)
    U=mc.model.U
    N_slices=mc.parameters.slices
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]

    ϕ0_bar= zeros(Float64, N_slices)
    for slice in 1:N_slices
        for i in 1:N
            iy, ix=fldmod1(i, L)                                    
            if iseven(ix+iy)
                Δpi=((-1)^(ix-1)+(-1)^(iy-1))
                temp0=zero(Float64)
                for k in 1:N
                    kPi=Bsrcdir2trg[k, i]           # k + i  
                    temp0 += sum(ϕ_field[: ,kPi ,slice] .*ϕ_field[:, k,slice])
                end
                ϕ0_bar[slice] +=Δpi * temp0
            end
        end
    end
    return ϕ0_bar ./ (2U*N^2)
end

"""
Calculates ̄Φ^(A₁')_ϕ (ℓ), returning the vector ℓ=1,…,N_τ
"""
function calc_ΦA1p_bar(mc ::DQMC, ϕ_field::Array)
    Nϕ=size(ϕ_field, 1)
    U=mc.model.U
    N_slices=mc.parameters.slices
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]

    ϕA1p_bar= zeros(Float64, N_slices)
    for slice in 1:N_slices
        for i in 1:N
            iy, ix=fldmod1(i, L)                                    
            if iseven(ix+iy)
                Δpi=((-1)^(ix-1)+(-1)^(iy-1))
                temp0=zero(Float64)
                for k in 1:N
                    ky, kx=fldmod1(k, L)   
                    kPi=Bsrcdir2trg[k, i]           # k + i  
                    temp0 += (-1)^(kx+ky) *sum(ϕ_field[: ,kPi ,slice] .*ϕ_field[:, k,slice])
                end
                ϕA1p_bar[slice] +=Δpi * temp0
            end
        end
    end
    return ϕA1p_bar ./ (2U*N^2)
end

"""
Calculates ̄Φ^(B₁')_ϕ (ℓ), returning the matrix, η=1,…,Nϕ, ℓ=1,…,N_τ
"""
function calc_ΦB1p_bar(mc ::DQMC, ϕ_field::Array)
    Nϕ=size(ϕ_field, 1)
    U=mc.model.U
    N_slices=mc.parameters.slices
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]
    if Nϕ==1
        println("Φ^(B₁') is not defined for the Ising case with Nϕ=1")
        return nothing
    end
    if Nϕ==2
        ϕB1p_bar= zeros(Float64, 1, N_slices)
        for slice in 1:N_slices
            for i in 1:N
                iy, ix=fldmod1(i, L)                                    
                if isodd(ix+iy)
                    Δmi=((-1)^(ix-1)-(-1)^(iy-1))
                    temp0=zero(Float64)
                    for k in 1:N
                        ky, kx=fldmod1(k, L)   
                        kPi=Bsrcdir2trg[k, i]           # k + i  
                        temp0 += (-1)^(kx+ky) *(ϕ_field[1 ,kPi ,slice] *ϕ_field[2, k,slice]-ϕ_field[2 ,kPi ,slice] *ϕ_field[1, k,slice])
                    end
                    ϕB1p_bar[1, slice] +=Δmi * temp0
                end
            end
        end
        return ϕB1p_bar ./ (2U*N^2)
    else
        ϕB1p_bar= zeros(Float64, 3, N_slices)
        for slice in 1:N_slices
            for i in 1:N
                iy, ix=fldmod1(i, L)                                    
                if isodd(ix+iy)
                    Δmi=((-1)^(ix-1)-(-1)^(iy-1))
                    temp0=zeros(Float64, 3)
                    for k in 1:N
                        ky, kx=fldmod1(k, L)   
                        kPi=Bsrcdir2trg[k, i]           # k + i  
                        temp0[1] += (-1)^(kx+ky) *(ϕ_field[2 ,kPi ,slice] *ϕ_field[3, k,slice]-ϕ_field[3 ,kPi ,slice] *ϕ_field[2, k,slice])
                        temp0[2] += (-1)^(kx+ky) *(ϕ_field[3 ,kPi ,slice] *ϕ_field[1, k,slice]-ϕ_field[1 ,kPi ,slice] *ϕ_field[3, k,slice])
                        temp0[3] += (-1)^(kx+ky) *(ϕ_field[1 ,kPi ,slice] *ϕ_field[2, k,slice]-ϕ_field[2 ,kPi ,slice] *ϕ_field[1, k,slice])
                    end
                    ϕB1p_bar[:, slice] .+=Δmi .* temp0[:]
                end
            end
        end
        return ϕB1p_bar ./ (2U*N^2)
    end
end