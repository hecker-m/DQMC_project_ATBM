
abstract type AbstractMagnBosonField <: AbstractField end
abstract type AbstractContMBF <: AbstractMagnBosonField end

# struct MagnBosonField <: AbstractContMBF
#     temp_conf::Array{Float64}
#     conf::Array{Float64}
#     delta_tau_copy::Float64     #for calculation of energy_boson it proves useful to have it saved here too
#     temp_vec::Vector{Float64}
# end

# function MagnBosonField(param::DQMCParameters,model::Model, Nϕ::Int64= param.Nϕ)
#     MagnBosonField(Array{Float64}(undef, Nϕ, length(lattice(model)), param.slices),
#     Array{Float64}(undef, Nϕ,length(lattice(model)), param.slices), 
#     param.delta_tau, zeros(Float64,Nϕ))
# end

struct Cont_MBF1 <: AbstractContMBF #Continuous magnetic boson field with Nϕ=1
    temp_conf::Array{Float64}
    conf::Array{Float64}
    delta_tau_copy::Float64     #for calculation of energy_boson it proves useful to have it saved here too
    α::Float64                  #α=  sqrt(2U) * δτ, useful for interaction matrix
    temp_vec::Vector{Float64}
end
function Cont_MBF1(param::DQMCParameters,model::Model)

    Cont_MBF1(Array{Float64}(undef, 1, length(lattice(model)), param.slices),
    Array{Float64}(undef, 1,length(lattice(model)), param.slices), 
    param.delta_tau, sqrt(2model.U) * param.delta_tau, Vector{Float64}(undef,1))
end
struct Cont_MBF2 <: AbstractContMBF     #Continuous magnetic boson field with Nϕ=2
    temp_conf::Array{Float64}
    conf::Array{Float64}
    delta_tau_copy::Float64     #for calculation of energy_boson it proves useful to have it saved here too
    α::Float64                  #α=  sqrt(2U) * δτ, useful for interaction matrix
    temp_vec::Vector{Float64}
end
function Cont_MBF2(param::DQMCParameters,model::Model)
    Cont_MBF2(Array{Float64}(undef, 2, length(lattice(model)), param.slices),
    Array{Float64}(undef, 2,length(lattice(model)), param.slices), 
    param.delta_tau, sqrt(2model.U) * param.delta_tau, Vector{Float64}(undef,2))
end
struct Cont_MBF3 <: AbstractContMBF     #Continuous magnetic boson field with Nϕ=3
    temp_conf::Array{Float64}
    conf::Array{Float64}
    delta_tau_copy::Float64     #for calculation of energy_boson it proves useful to have it saved here too
    α::Float64                  # α= sqrt(2U) * δτ, useful for interaction matrix
    temp_vec::Vector{Float64}
end
function Cont_MBF3(param::DQMCParameters,model::Model)

    Cont_MBF3(Array{Float64}(undef, 3, length(lattice(model)), param.slices),
    Array{Float64}(undef, 3,length(lattice(model)), param.slices), 
    param.delta_tau, sqrt(2model.U) * param.delta_tau, Vector{Float64}(undef,3))
end
unique_flavors(::AbstractContMBF) = 4
#energy_boson(::MagnBosonField, conf=nothing) = 0.0
energy_boson(f::AbstractContMBF, conf = f.conf) =0.5*f.delta_tau_copy*sum(conf .^2)


interaction_eltype(f::Cont_MBF1) =Float64 
interaction_eltype(f::Cont_MBF2) =ComplexF64 
interaction_eltype(f::Cont_MBF3) =ComplexF64 
#interaction_matrix_type(f::AbstractMagnBosonField, ::Model) = SparseMatrixCSC{interaction_eltype(f), Int64}
interaction_matrix_type(f::AbstractMagnBosonField, m::Model) = Matrix{m.peierls ==true ? ComplexF64 : interaction_eltype(f)}

function init_interaction_matrix(f::AbstractMagnBosonField, m::Model)
    flv = unique_flavors(f, m)
    N = length(lattice(m))
    #spzeros(interaction_eltype(f),N*flv, N*flv)
    zeros(m.peierls ==true ? ComplexF64 : interaction_eltype(f),N*flv, N*flv)
end

vmul!(C::Matrix, A::Hermitian, B::SparseMatrixCSC) = vmul!(C, A.data, B) 
vmul!(C::Matrix, A::SparseMatrixCSC, B::Hermitian) = vmul!(C, A, B.data) 
function vmul!(C::Matrix, A::SparseMatrixCSC, B::Matrix)
    @debug "vmul!($(typeof(C)), $(typeof(A)), $(typeof(B)))"
    mul!(C, A, B)        
end
function vmul!(C::Matrix, A::Matrix, B::SparseMatrixCSC)
    @debug "vmul!($(typeof(C)), $(typeof(A)), $(typeof(B)))"
    mul!(C, A, B)        
end

Base.rand(f::AbstractContMBF) = rand(size(f.conf)[1], size(f.conf)[2],size(f.conf)[3])
function Random.rand!(f::AbstractContMBF; range::Float64=1.0)
    rand!(f.conf)
    f.conf .= 2range .* f.conf .-range  #intial configuration is ∈[-range,range]
end
compressed_conf_type(::AbstractContMBF) = Array
compressed_conf_type(::Type{<: DQMC}, ::Type{<: TwoBandModel}) = Array

function ConfigRecorder(::Type{<: AbstractMagnBosonField}, rate::Int = 10)
    ConfigRecorder{Array}(rate)
end
function ConfigRecorder(::Type{<: AbstractContMBF}, rate::Int = 10)
    ConfigRecorder{Array}(rate)
end

function Base.push!(c::ConfigRecorder, field::AbstractContMBF, sweep)
    (sweep % c.rate == 0) && push!(c.configs, copy(field.conf)) 
    nothing
end

function BufferedConfigRecorder(::Type{<: AbstractContMBF}, filename; rate = 10, chunk_size = 1000)
    BufferedConfigRecorder{Array}(filename, rate, chunk_size)
end

function Base.push!(cr::BufferedConfigRecorder, field::AbstractContMBF, sweep)
    (sweep % cr.rate == 0) && _push!(cr, field.conf)
    nothing
end

# cosmetics
import Base.summary
import Base.show
function Base.summary(field::AbstractContMBF)
    " Continuous magnetic bosonic field"
end
function Base.show(io::IO, field::AbstractContMBF)
    println(io, " Continuous magnetic bosonic field")
    println(io, "\tNumber of components Nϕ= $(size(field.conf)[1])")
end


"""
Draw a random number from a uniform distribution over an interval [-b/2, b/2].
"""
@inline randuniform(b::Float64) = -b/2 + b * rand() # taken from Distributions.jl: https://tinyurl.com/ycr9jnt4
randuniform(b::Float64, d::Int) = begin
    x = Vector{Float64}(undef, d)
    @inbounds for k in Base.OneTo(d)
      x[k] = randuniform(b)
    end
    x
end




@inline function calculate_detratio!(mc::DQMC, model::Model, i::Int, 
    new_op::Vector{G}) where {G<:Number}
    s=mc.stack
    N=length(model.l)
    interaction_matrix_exp_op!(mc, mc.model, mc.field.conf[:,i,s.current_slice], -1., s.field_cache.eVop1) #V1i
    interaction_matrix_exp_op!(mc, mc.model, new_op, 1.,s.field_cache.eVop2) #V2i

    mul!(s.field_cache.eVop3, s.field_cache.eVop1, s.field_cache.eVop2)
    s.field_cache.Δ .= s.field_cache.eVop3 .- s.eye_flv

    s.field_cache.eVop3 .= s.eye_flv .- s.greens[i:N:end,i:N:end]
    mul!(s.field_cache.eVop1, s.field_cache.eVop3, s.field_cache.Δ)
    s.field_cache.R .= s.eye_flv .+ s.field_cache.eVop1 # K_tilde=I+(I-G)*Δ is stored in R
    return det(s.field_cache.R)
end

@inline function propose_local(mc::DQMC, model::Model, f::Cont_MBF1, i::Int, slice::Int)
    f.temp_vec[1]=f.conf[1,i, slice] + randuniform(mc.parameters.box_local)   
    ΔE_boson ::Float64=0.5*mc.parameters.delta_tau*(f.temp_vec[1]^2 - f.conf[1,i, slice]^2)
    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return detratio, ΔE_boson, f.temp_vec
end

@inline function propose_local(mc::DQMC, model::Model, f::AbstractContMBF, i::Int, slice::Int)
    f.temp_vec[:] = f.conf[:,i, slice] + randuniform(mc.parameters.box_local, mc.parameters.Nϕ)
    ΔE_boson ::Float64=0.5*mc.parameters.delta_tau*(f.temp_vec⋅f.temp_vec-f.conf[:,i, slice]⋅f.conf[:,i, slice])
    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return detratio, ΔE_boson, f.temp_vec
end

# accept local works for continuous and discrete fields alike
@inline @bm function accept_local!(mc::DQMC, m::Model, f::AbstractMagnBosonField,
     i::Int, slice::Int,  detratio, ΔE_boson, new_op)
    update_greens!(mc.stack.field_cache, mc.stack.greens, i, length(mc.model.l))
    @inbounds f.conf[:, i, slice] = new_op
    return nothing
end







#= already defined in field.jl
function update_greens!(cache::StandardFieldCache, G, i, N)
    # calculate Δ R⁻¹
    vldiv22!(cache, cache.R, cache.Δ)
    
    # calculate (I - G)[:, i:N:end]
    vsub!(cache.IG, I, G, i, N)

    # calculate {Δ R⁻¹} * G[i:N:end, :]
    vmul!(cache.G, cache.invRΔ, G, i, N)

    # update greens function 
    # G[m, n] -= {(I - G)[m, i:N:end]} {{Δ R⁻¹} * G[i:N:end, n]}
    vsubkron!(G, cache.IG, cache.G)

    nothing
end
=#


###############################################
## interactin matrix exponentials (continuous fields)
##############################################


# exp(-power*δτ*V)
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, f::Cont_MBF1, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N

        ch3 = cosh(f.α *f.conf[1,i,slice])  # α= sqrt(2U) *δτ
        sh3 = sinh(f.α *f.conf[1,i,slice])
   
        result[i,i]         = ch3
        result[i+N,i+N]     = ch3
        result[i+2N,i+2N]   = ch3
        result[i+3N,i+3N]   = ch3

        result[i,i+2N]      = power*sh3
        result[i+N,i+3N]    = -power*sh3
        result[i+2N,i]      = power*sh3
        result[i+3N,i+N]    = -power*sh3
    end
    nothing
end

@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, f::Cont_MBF2, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N

        ch1 = cosh(f.α *f.conf[1,i,slice])   # α= sqrt(2U) *δτ
        sh1 = sinh(f.α *f.conf[1,i,slice])
        ch2 = cosh(f.α *f.conf[2,i,slice])
        sh2 = sinh(f.α *f.conf[2,i,slice])
   
        result[i,i]         = ch1*ch2 + im*power*sh1*sh2
        result[i+N,i+N]     = ch1*ch2 - im*power*sh1*sh2
        result[i+2N,i+2N]   = ch1*ch2 + im*power*sh1*sh2
        result[i+3N,i+3N]   = ch1*ch2 - im*power*sh1*sh2

        result[i,i+3N]      = power*(ch2*sh1 - im*ch1*sh2)
        result[i+N,i+2N]    = power*(ch2*sh1 + im*ch1*sh2)
        result[i+2N,i+N]    = power*(ch2*sh1 - im*ch1*sh2)
        result[i+3N,i]      = power*(ch2*sh1 + im*ch1*sh2)

    end
    nothing
end

#exact function
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, f::Cont_MBF3, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N

        ch1 = cosh(f.α *f.conf[1,i,slice])   # α= sqrt(2U) *δτ
        sh1 = sinh(f.α *f.conf[1,i,slice])
        ch2 = cosh(f.α *f.conf[2,i,slice])
        sh2 = sinh(f.α *f.conf[2,i,slice])
        ch3 = cosh(f.α *f.conf[3,i,slice])
        sh3 = sinh(f.α *f.conf[3,i,slice])
   
        result[i,i]         = ch1*ch2*ch3 + im*power*sh1*sh2*ch3
        result[i+N,i+N]     = ch1*ch2*ch3 - im*power*sh1*sh2*ch3
        result[i+2N,i+2N]   = ch1*ch2*ch3 + im*power*sh1*sh2*ch3
        result[i+3N,i+3N]   = ch1*ch2*ch3 - im*power*sh1*sh2*ch3

        result[i,i+N]       = power*(-ch2*sh1*sh3+im*ch1*sh2*sh3)
        result[i+N,i]       = power*(ch2*sh1*sh3+im*ch1*sh2*sh3)
        result[i+2N,i+3N]   = power*(-ch2*sh1*sh3+im*ch1*sh2*sh3)
        result[i+3N,i+2N]   = power*(ch2*sh1*sh3+im*ch1*sh2*sh3)

        result[i,i+2N]      = power*ch1*ch2*sh3+ im*sh1*sh2*sh3
        result[i+N,i+3N]    = -power*ch1*ch2*sh3+ im*sh1*sh2*sh3
        result[i+2N,i]      = power*ch1*ch2*sh3+ im*sh1*sh2*sh3
        result[i+3N,i+N]    = -power*ch1*ch2*sh3+ im*sh1*sh2*sh3

        result[i,i+3N]      = power*ch3*(ch2*sh1 - im*ch1*sh2)
        result[i+N,i+2N]    = power*ch3*(ch2*sh1 + im*ch1*sh2)
        result[i+2N,i+N]    = power*ch3*(ch2*sh1 - im*ch1*sh2)
        result[i+3N,i]      = power*ch3*(ch2*sh1 + im*ch1*sh2)
    end
    return nothing
end

# exp(-power*δτ*V)
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel, op::Vector{T}, 
    power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where {T<:Number,G}
    interaction_matrix_exp_op!(mc, model, field(mc), op, power, eVop)
    return nothing
end
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel, f::Cont_MBF1, op::Vector{Float64}, 
    power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    ch3 = cosh(f.α * op[1]) # α= sqrt(2U) *δτ
    sh3 = sinh(f.α * op[1])

    eVop[1,1] = ch3
    eVop[2,2] = ch3
    eVop[3,3] = ch3
    eVop[4,4] = ch3
    eVop[1,2] = zero(G)
    eVop[2,1] = zero(G)
    eVop[3,4] = zero(G)
    eVop[4,3] = zero(G)
    eVop[1,3] = power*sh3
    eVop[2,4] = -power*sh3
    eVop[3,1] = power*sh3
    eVop[4,2] = -power*sh3
    eVop[1,4] = zero(G)
    eVop[2,3] = zero(G)
    eVop[3,2] = zero(G)
    eVop[4,1] = zero(G)
    return nothing
end

@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel,  f::Cont_MBF2, op::Vector{Float64}, 
    power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    ch1 = cosh(f.α * op[1])  # α= sqrt(2U) *δτ
    sh1 = sinh(f.α * op[1])
    ch2 = cosh(f.α * op[2])
    sh2 = sinh(f.α * op[2])

    eVop[1,1] = ch1*ch2+im*power*sh1*sh2
    eVop[2,2] = ch1*ch2-im*power*sh1*sh2
    eVop[3,3] = ch1*ch2+im*power*sh1*sh2
    eVop[4,4] = ch1*ch2-im*power*sh1*sh2
    eVop[1,2] = zero(G)
    eVop[2,1] = zero(G)
    eVop[3,4] = zero(G)
    eVop[4,3] = zero(G)
    eVop[1,3] = zero(G)
    eVop[2,4] = zero(G)
    eVop[3,1] = zero(G)
    eVop[4,2] = zero(G)
    eVop[1,4] = power*(ch2*sh1 - im*ch1*sh2)
    eVop[2,3] = power*(ch2*sh1 + im*ch1*sh2)
    eVop[3,2] = power*(ch2*sh1 - im*ch1*sh2)
    eVop[4,1] = power*(ch2*sh1 + im*ch1*sh2)
    return nothing     
end

#exact function 
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel,  f::Cont_MBF3, op::Vector{Float64},
    power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    ch1 = cosh(f.α * op[1])  # α= sqrt(2U) *δτ
    sh1 = sinh(f.α * op[1])
    ch2 = cosh(f.α * op[2])
    sh2 = sinh(f.α * op[2])
    ch3 = cosh(f.α * op[3])
    sh3 = sinh(f.α * op[3])

    eVop[1,1] = ch1*ch2*ch3+im*power*sh1*sh2*ch3
    eVop[2,2] = ch1*ch2*ch3-im*power*sh1*sh2*ch3
    eVop[3,3] = ch1*ch2*ch3+im*power*sh1*sh2*ch3
    eVop[4,4] = ch1*ch2*ch3-im*power*sh1*sh2*ch3
    eVop[1,2] = power*(-ch2*sh1*sh3+im*ch1*sh2*sh3)
    eVop[2,1] = power*(ch2*sh1*sh3+im*ch1*sh2*sh3)
    eVop[3,4] = power*(-ch2*sh1*sh3+im*ch1*sh2*sh3)
    eVop[4,3] = power*(ch2*sh1*sh3+im*ch1*sh2*sh3)
    eVop[1,3] = power*ch1*ch2*sh3+ im*sh1*sh2*sh3
    eVop[2,4] = -power*ch1*ch2*sh3+ im*sh1*sh2*sh3
    eVop[3,1] = power*ch1*ch2*sh3+ im*sh1*sh2*sh3
    eVop[4,2] = -power*ch1*ch2*sh3+ im*sh1*sh2*sh3
    eVop[1,4] = power*ch3*(ch2*sh1 - im*ch1*sh2)
    eVop[2,3] = power*ch3*(ch2*sh1 + im*ch1*sh2)
    eVop[3,2] = power*ch3*(ch2*sh1 - im*ch1*sh2)
    eVop[4,1] = power*ch3*(ch2*sh1 + im*ch1*sh2)
    return nothing
end

