
abstract type AbstractDiscreteMBF <: AbstractMagnBosonField end


"""
Discrete_MBF1

    IsingZ-magnetic boson field with Nϕ=1

    Note that this is the full version, i.e. it is not symmetry-optimized.
"""
struct Discrete_MBF1 <: AbstractDiscreteMBF #Continuous magnetic boson field with Nϕ=1
    α::Float64          # α=sqrt(δτ*U)
    w::Vector{Float64}  #holds the four bosonic weights
    x::Vector{Float64}  #holds the four exponential weights
    choices::Matrix{Int8}
    temp_conf::Array{Int8, 3} #stores numbers η ∈ {1,2,3,4}
    conf::Array{Int8, 3}
    temp_vec::Vector{Int8}
end
function Discrete_MBF1(param::DQMCParameters, model::Model)
    α = sqrt(param.delta_tau *model.U)
    s6 = sqrt(6)
    ws = Float64[1 - s6/3, 1 + s6/3, 1 + s6/3, 1 - s6/3]
    xs = Float64[-sqrt(6 + 2s6), -sqrt(6 - 2s6), sqrt(6 - 2s6), sqrt(6 + 2s6)]
    choices = Int8[2 3 4; 1 3 4; 1 2 4; 1 2 3]
    Discrete_MBF1( α, ws, xs, choices,
        Array{Int8}(undef, 1, length(lattice(model)), param.slices),
        Array{Int8}(undef, 1, length(lattice(model)), param.slices), 
        Vector{Int8}(undef,1))
end
"""
Discrete_MBF1_symm

    Symmetry-optimized IsingZ-magnetic boson field with Nϕ=1
"""
struct Discrete_MBF1_symm <: AbstractDiscreteMBF #Continuous IsingZ-magnetic boson field with Nϕ=1
    α::Float64          # α=sqrt(δτ*U)
    w::Vector{Float64}  #holds the four bosonic weights
    x::Vector{Float64}  #holds the four exponential weights
    choices::Matrix{Int8}
    temp_conf::Array{Int8, 3} #stores numbers η ∈ {1,2,3,4}
    conf::Array{Int8, 3}
    temp_vec::Vector{Int8}
end

function Discrete_MBF1_symm(param::DQMCParameters, model::Model)
    α = sqrt(param.delta_tau *model.U)
    s6 = sqrt(6)
    ws = Float64[1 - s6/3, 1 + s6/3, 1 + s6/3, 1 - s6/3]
    xs = Float64[-sqrt(6 + 2s6), -sqrt(6 - 2s6), sqrt(6 - 2s6), sqrt(6 + 2s6)]
    choices = Int8[2 3 4; 1 3 4; 1 2 4; 1 2 3]
    Discrete_MBF1_symm( α, ws, xs, choices,
        Array{Int8}(undef, 1, length(lattice(model)), param.slices),
        Array{Int8}(undef, 1, length(lattice(model)), param.slices), 
        Vector{Int8}(undef,1))
end

"""
Discrete_MBF1_X

    Ising>X-magnetic boson field with Nϕ=1

    Note that this is the full version, i.e. it is not symmetry-optimized.
"""
struct Discrete_MBF1_X <: AbstractDiscreteMBF #Continuous magnetic boson field with Nϕ=1
    α::Float64          # α=sqrt(δτ*U)
    w::Vector{Float64}  #holds the four bosonic weights
    x::Vector{Float64}  #holds the four exponential weights
    choices::Matrix{Int8}
    temp_conf::Array{Int8, 3} #stores numbers η ∈ {1,2,3,4}
    conf::Array{Int8, 3}
    temp_vec::Vector{Int8}
end
function Discrete_MBF1_X(param::DQMCParameters, model::Model)
    α = sqrt(param.delta_tau *model.U)
    s6 = sqrt(6)
    ws = Float64[1 - s6/3, 1 + s6/3, 1 + s6/3, 1 - s6/3]
    xs = Float64[-sqrt(6 + 2s6), -sqrt(6 - 2s6), sqrt(6 - 2s6), sqrt(6 + 2s6)]
    choices = Int8[2 3 4; 1 3 4; 1 2 4; 1 2 3]
    Discrete_MBF1_X( α, ws, xs, choices,
        Array{Int8}(undef, 1, length(lattice(model)), param.slices),
        Array{Int8}(undef, 1, length(lattice(model)), param.slices), 
        Vector{Int8}(undef,1))
end

"""
Discrete_MBF1_X_symm

    Symmetry-optimized IsingX-magnetic boson field with Nϕ=1
"""
struct Discrete_MBF1_X_symm <: AbstractDiscreteMBF #Continuous magnetic boson field with Nϕ=1
    α::Float64          # α=sqrt(δτ*U)
    w::Vector{Float64}  #holds the four bosonic weights
    x::Vector{Float64}  #holds the four exponential weights
    choices::Matrix{Int8}
    temp_conf::Array{Int8, 3} #stores numbers η ∈ {1,2,3,4}
    conf::Array{Int8, 3}
    temp_vec::Vector{Int8}
end
function Discrete_MBF1_X_symm(param::DQMCParameters, model::Model)
    α = sqrt(param.delta_tau *model.U)
    s6 = sqrt(6)
    ws = Float64[1 - s6/3, 1 + s6/3, 1 + s6/3, 1 - s6/3]
    xs = Float64[-sqrt(6 + 2s6), -sqrt(6 - 2s6), sqrt(6 - 2s6), sqrt(6 + 2s6)]
    choices = Int8[2 3 4; 1 3 4; 1 2 4; 1 2 3]
    Discrete_MBF1_X_symm( α, ws, xs, choices,
        Array{Int8}(undef, 1, length(lattice(model)), param.slices),
        Array{Int8}(undef, 1, length(lattice(model)), param.slices), 
        Vector{Int8}(undef,1))
end

"""
Discrete_MBF2

    XY-magnetic boson field with Nϕ=2
    Note that this is the full version, i.e. it is not symmetry-optimized.
"""
struct Discrete_MBF2 <: AbstractDiscreteMBF #Continuous magnetic boson field with Nϕ=2
    α::Float64
    w::Vector{Float64}
    x::Vector{Float64}
    choices::Matrix{Vector{Int8}}
    temp_conf::Array{Int8, 3} #stores numbers η ∈ {1,2,3,4}
    conf::Array{Int8, 3}
    temp_vec::Vector{Int8}
end

function Discrete_MBF2(param::DQMCParameters, model::Model)
    α = sqrt(param.delta_tau *model.U)
    s6 = sqrt(6)
    ws = Float64[1 - s6/3, 1 + s6/3, 1 + s6/3, 1 - s6/3]
    xs = Float64[-sqrt(6 + 2s6), -sqrt(6 - 2s6), sqrt(6 - 2s6), sqrt(6 + 2s6)]
    choices = [[Int8(i), Int8(j)] for i ∈ 1:4, j ∈ 1:4]
    Discrete_MBF2( α, ws, xs, choices,
        Array{Int8}(undef, 2, length(lattice(model)), param.slices),
        Array{Int8}(undef, 2, length(lattice(model)), param.slices), 
        Vector{Int8}(undef,2))
end
"""
Discrete_MBF3

    Heisenberg-magnetic boson field with Nϕ=3
"""
struct Discrete_MBF3 <: AbstractDiscreteMBF #Continuous magnetic boson field with Nϕ=3
    α::Float64
    w::Vector{Float64}
    x::Vector{Float64}
    #choices::Array{Vector{Int8}, 3}
    choices::Array{Int8, 5}
    temp_conf::Array{Int8, 3} #stores numbers η ∈ {1,2,3,4}
    conf::Array{Int8, 3}
    temp_vec::Vector{Int8}
end

function make_choices(xs)
    fac0=1;
    mat2=Array{Int8, 5}(undef, 3, 63+31*(fac0-1) ,4, 4, 4)
    for i0 in 1:4, j0 in 1:4, k0 in 1:4
        μ0=[i0, j0, k0];
        count=0;
        for i in 1:4, j in 1:4, k in 1:4
            vec=[i, j, k]
            if vec != μ0
                xs[μ0] ⋅ xs[vec]>0 ? fac=fac0 : fac=1
                for f in 1:fac
                    count+=1;
                    mat2[:, count, i0, j0, k0]=vec
                end
            end
        end
    end
    return mat2
end

function Discrete_MBF3(param::DQMCParameters, model::Model)
    α = sqrt(param.delta_tau *model.U)
    s6 = sqrt(6)
    ws = Float64[1 - s6/3, 1 + s6/3, 1 + s6/3, 1 - s6/3]
    xs = Float64[-sqrt(6 + 2s6), -sqrt(6 - 2s6), sqrt(6 - 2s6), sqrt(6 + 2s6)]
    # choices = [[Int8(i), Int8(j), Int8(k)] for i ∈ 1:4, j ∈ 1:4, k ∈ 1:4]
    choices = make_choices(xs)
    Discrete_MBF3( α, ws, xs, choices,
        Array{Int8}(undef, 3, length(lattice(model)), param.slices),
        Array{Int8}(undef, 3, length(lattice(model)), param.slices),
        Vector{Int8}(undef,3))
end



unique_flavors(::AbstractDiscreteMBF) = 4

interaction_eltype(f::Union{Discrete_MBF1,Discrete_MBF1_X}) =Float64 
interaction_eltype(f::Discrete_MBF2) =ComplexF64 
interaction_eltype(f::Discrete_MBF3) =ComplexF64 


# These represent (-2, -1, +1, 2)
const _DMBF_VALS = (Int8(1), Int8(2), Int8(3), Int8(4))
Base.rand(f::AbstractDiscreteMBF) = rand(_DMBF_VALS, size(f.conf))
Random.rand!(f::AbstractDiscreteMBF) = rand!(f.conf, _DMBF_VALS)
Random.rand!(f::AbstractMagnBosonField; kwargs...)=rand!(f)


################################################################################
### Utility
################################################################################

compressed_conf_type(::AbstractDiscreteMBF) = BitArray
function compress(f::AbstractDiscreteMBF)
    # converts (1, 2, 3, 4) -> (00, 01, 10, 11)
    BitArray((div(v-1, 2), (v-1) % 2)[step] for v in f.conf for step in (1, 2))
end
function decompress(f::AbstractDiscreteMBF, c)
    # converts (00, 01, 10, 11) -> (1, 2, 3, 4)
    conf = similar(f.conf)
    for i in eachindex(conf)
        #  1    +    2 * bit1    +    bit2
        conf[i] = Int8(1) + Int8(2) * c[2i-1] + Int8(c[2i])
    end
    conf
end
function decompress!(f::AbstractDiscreteMBF, c)
    for i in eachindex(f.conf)
        #  1    +    2 * bit1    +    bit2
        f.conf[i]=Int8(1) + Int8(2) * c[i] + Int8(c[i+1]) 
        #I added f.conf[i]=
        #Not sure if it worked before
    end
    #f.conf
end

function ConfigRecorder(::Type{<: AbstractDiscreteMBF}, rate::Int = 10)
    ConfigRecorder{BitArray}(rate)
end

function Base.push!(c::ConfigRecorder, field::AbstractDiscreteMBF, sweep)
    (sweep % c.rate == 0) && push!(c.configs, compress(field))
    nothing
end

function BufferedConfigRecorder(::Type{<: AbstractDiscreteMBF}, 
        filename; rate = 10, chunk_size = 1000)
    BufferedConfigRecorder{BitArray}(filename, rate, chunk_size)
end

function Base.push!(cr::BufferedConfigRecorder, field::AbstractDiscreteMBF, sweep)
    (sweep % cr.rate == 0) && _push!(cr, compress(field))
    nothing
end


# cosmetics
import Base.summary
import Base.show
function Base.summary(field::AbstractDiscreteMBF)
    " Discrete magnetic bosonic field"
end
function Base.show(io::IO, field::AbstractDiscreteMBF)
    println(io, " Discrete magnetic bosonic field")
    println(io, "\tNumber of components Nϕ= $(size(field.conf)[1])")
end




@inline function propose_local(mc::DQMC, model::Model, 
    f::Union{Discrete_MBF1, Discrete_MBF1_X}, i::Int, slice::Int)
    η_old = f.conf[1,i, slice]

    f.temp_vec[:]=@inbounds [f.choices[η_old, rand(1:3)]]

    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return detratio * f.w[f.temp_vec[1]]/f.w[η_old], 0.0, f.temp_vec
end

@inline function propose_local(mc::DQMC, model::Model, f::Discrete_MBF2, i::Int, slice::Int)

    rand1=rand(1:15)
    px=mod1(f.conf[1,i, slice] + mod(rand1, 4),4)
    py=mod1(f.conf[2,i, slice] + div(rand1, 4),4)
    f.temp_vec[:]=@inbounds f.choices[px, py]

    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return detratio * (f.w[f.temp_vec[1]]*f.w[f.temp_vec[2]])/
        (f.w[f.conf[1,i, slice]]*f.w[f.conf[2,i, slice]]), 
            0.0, f.temp_vec
end

@inline function propose_local(mc::DQMC, model::Model, f::Discrete_MBF3, i::Int, slice::Int)
    # rand1=rand(1:63)
    # px=mod1(f.conf[1,i, slice] + mod(rand1, 4), 4)
    # py=mod1(f.conf[2,i, slice] + div(mod(rand1, 16), 4), 4)
    # pz=mod1(f.conf[3,i, slice] + div(rand1, 16), 4)
    # f.temp_vec[:]=@inbounds f.choices[px, py, pz]
    rand1=rand(1:size(f.choices)[2])
    ix=f.conf[1,i, slice]
    iy=f.conf[2,i, slice]
    iz=f.conf[3,i, slice]
    f.temp_vec[:]=@inbounds f.choices[:, rand1, ix, iy, iz]
    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return detratio * (f.w[f.temp_vec[1]]*f.w[f.temp_vec[2]]*f.w[f.temp_vec[3]])/
        (f.w[f.conf[1,i, slice]]*f.w[f.conf[2,i, slice]]*f.w[f.conf[3,i, slice]]), 
            0.0, f.temp_vec
end
####
# fallback function for the Discrete MBF global updates
####
@bm function propose_global_from_conf(mc::DQMC, m::Model, f::AbstractDiscreteMBF)
    η_new=conf(f)
    η_old=temp_conf(f)

    @assert mc.stack.current_slice == 1     
    @assert mc.stack.direction == 1
    mc.stack.tempvf .= mc.stack.Dl      
    inv_det(mc, current_slice(mc)-1, f)     #writes the (inverted) singular values into mc.stack.Dr, probably sorted ?

    detratio = 1.0
    for i in eachindex(mc.stack.tempvf)
        detratio *= mc.stack.tempvf[i] * mc.stack.Dr[i]
    end

    fac=1.0;
    for i in eachindex(η_new) 
        fac *= f.w[η_new[i]]/f.w[η_old[i]]
    end
    return detratio * fac, 0.0, nothing
end

###############################################
## interactin matrix exponentials (discrete fields)
##############################################

# exp(-power*δτ*V)
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel, f::Discrete_MBF1,
    η_vec::Vector{Int8}, power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    sh3 = sinh(f.α * f.x[η_vec[1]]) # α=sqrt(δτ*U)
    ch3 = cosh(f.α * f.x[η_vec[1]])

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
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel, f::Discrete_MBF1_X, 
    η_vec::Vector{Int8}, power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    sh1 = sinh(f.α * f.x[η_vec[1]])   # α=sqrt(δτ*U)
    ch1 = cosh(f.α * f.x[η_vec[1]])

    eVop[1,1] = ch1
    eVop[2,2] = ch1
    eVop[3,3] = ch1
    eVop[4,4] = ch1
    eVop[1,2] = zero(G)
    eVop[2,1] = zero(G)
    eVop[3,4] = zero(G)
    eVop[4,3] = zero(G)
    eVop[1,3] = zero(G)
    eVop[2,4] = zero(G)
    eVop[3,1] = zero(G)
    eVop[4,2] = zero(G)
    eVop[1,4] = power*sh1
    eVop[2,3] = power*sh1
    eVop[3,2] = power*sh1
    eVop[4,1] = power*sh1
    return nothing
end
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel, f::Discrete_MBF2, 
    η_vec::Vector{Int8}, power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    sh1 = sinh(f.α * f.x[η_vec[1]])   # α=sqrt(δτ*U)
    ch1 = cosh(f.α * f.x[η_vec[1]])
    sh2 = sinh(f.α * f.x[η_vec[2]])
    ch2 = cosh(f.α * f.x[η_vec[2]])

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
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel,f::Discrete_MBF3, 
    η_vec::Vector{Int8}, power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    sh1 = sinh(f.α * f.x[η_vec[1]])   # α=sqrt(δτ*U)
    ch1 = cosh(f.α * f.x[η_vec[1]])
    sh2 = sinh(f.α * f.x[η_vec[2]])
    ch2 = cosh(f.α * f.x[η_vec[2]])
    sh3 = sinh(f.α * f.x[η_vec[3]])
    ch3 = cosh(f.α * f.x[η_vec[3]])

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

# exp(-power*δτ*V) full (4N×4N) matrix
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, f::Discrete_MBF1, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N
        sh3 = sinh(f.α * f.x[f.conf[1,i,slice]])    # α=sqrt(δτ*U)
        ch3 = cosh(f.α * f.x[f.conf[1,i,slice]])

        result[i,i]         = ch3
        result[i+N,i+N]     = ch3
        result[i+2N,i+2N]   = ch3
        result[i+3N,i+3N]   = ch3

        result[i,i+2N]      = power*sh3
        result[i+N,i+3N]    = -power*sh3
        result[i+2N,i]      = power*sh3
        result[i+3N,i+N]    = -power*sh3
    end
    return nothing
end
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, f::Discrete_MBF1_X, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N
        
        sh1 = sinh(f.α * f.x[f.conf[1,i,slice]])  # α=sqrt(δτ*U)
        ch1 = cosh(f.α * f.x[f.conf[1,i,slice]])

        result[i,i]         = ch1
        result[i+N,i+N]     = ch1
        result[i+2N,i+2N]   = ch1
        result[i+3N,i+3N]   = ch1

        result[i,i+3N]      = power*sh1
        result[i+N,i+2N]    = power*sh1
        result[i+2N,i+N]    = power*sh1
        result[i+3N,i]      = power*sh1
    end
    return nothing
end
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, f::Discrete_MBF2, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N
        
        sh1 = sinh(f.α * f.x[f.conf[1,i,slice]])  # α=sqrt(δτ*U)
        ch1 = cosh(f.α * f.x[f.conf[1,i,slice]])
        sh2 = sinh(f.α * f.x[f.conf[2,i,slice]])
        ch2 = cosh(f.α * f.x[f.conf[2,i,slice]])

        result[i,i]         = ch1*ch2 + im*power*sh1*sh2
        result[i+N,i+N]     = ch1*ch2 - im*power*sh1*sh2
        result[i+2N,i+2N]   = ch1*ch2 + im*power*sh1*sh2
        result[i+3N,i+3N]   = ch1*ch2 - im*power*sh1*sh2

        result[i,i+3N]      = power*(ch2*sh1 - im*ch1*sh2)
        result[i+N,i+2N]    = power*(ch2*sh1 + im*ch1*sh2)
        result[i+2N,i+N]    = power*(ch2*sh1 - im*ch1*sh2)
        result[i+3N,i]      = power*(ch2*sh1 + im*ch1*sh2)
    end
    return nothing
end
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, f::Discrete_MBF3, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N
        sh1 = sinh(f.α * f.x[f.conf[1,i,slice]])  # α=sqrt(δτ*U)
        ch1 = cosh(f.α * f.x[f.conf[1,i,slice]])
        sh2 = sinh(f.α * f.x[f.conf[2,i,slice]])
        ch2 = cosh(f.α * f.x[f.conf[2,i,slice]])
        sh3 = sinh(f.α * f.x[f.conf[3,i,slice]])
        ch3 = cosh(f.α * f.x[f.conf[3,i,slice]])

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