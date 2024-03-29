


unique_flavors(::Union{Discrete_MBF1_symm, Abstract_DiscreteMBF1_X_symm, Cont_MBF1_X_symm}) = 2
interaction_eltype(::Union{Discrete_MBF1_symm, Abstract_DiscreteMBF1_X_symm, Cont_MBF1_X_symm}) =Float64 

unique_flavors(::Discrete_MBF2_symm) = 2
interaction_eltype(::Discrete_MBF2_symm) =ComplexF64 

function pad_to_unique_flavors(f::Discrete_MBF1_symm, m::TwoBandModel, mat)
    N = length(lattice(m))
    flv = unique_flavors(f, m)
    if size(mat, 1) == N * flv
        return mat
    elseif size(mat, 1) == 2flv*N
        hopp_mat=zeros(eltype(mat), N*flv, N*flv);
        hopp_mat[1:N, 1:N] .= mat[1:N, 1:N]
        hopp_mat[N+1:N+N, N+1:N+N] .= mat[2N+1:2N+N, 2N+1:2N+N]

        return hopp_mat
    else
        error("Failed to expand size $(size(mat)) matrix to size ($N * $flv, $N * $flv) ")
    end
end
function pad_to_unique_flavors(f::Union{Abstract_DiscreteMBF1_X_symm, Discrete_MBF2_symm, Cont_MBF1_X_symm}, 
        m::TwoBandModel, mat)
    N = length(lattice(m))
    flv = unique_flavors(f, m)
    if size(mat, 1) == N * flv
        return mat
    elseif size(mat, 1) == 2flv*N
        hopp_mat=zeros(eltype(mat), N*flv, N*flv);
        hopp_mat[1:N, 1:N] .= mat[1:N, 1:N]
        hopp_mat[N+1:N+N, N+1:N+N] .= mat[3N+1:3N+N, 3N+1:3N+N]

        return hopp_mat
    else
        error("Failed to expand size $(size(mat)) matrix to size ($N * $flv, $N * $flv) ")
    end
end


@inline function propose_local(mc::DQMC, model::Model, 
    f::Union{Discrete_MBF1_symm, Abstract_DiscreteMBF1_X_symm}, i::Int, slice::Int)

    k=size(f.choices)[2]

    η_old = f.conf[1,i, slice]
    f.temp_vec[:]=@inbounds [f.choices[η_old, rand(1:k)]]

    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return abs2(detratio) * f.w[f.temp_vec[1]]/f.w[η_old], 0.0, f.temp_vec
end

@inline function propose_local(mc::DQMC, model::Model, 
    f::Union{Discrete_8_MBF1_X_symm, Discrete_12_MBF1_X_symm, Discrete_16_MBF1_X_symm}, i::Int, slice::Int)

    k=size(f.choices)[2]

    η_old = f.conf[1,i, slice]
    f.temp_vec[:]=@inbounds [f.choices[η_old, rand(1:k)]]

    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return abs2(detratio) * f.w[f.temp_vec[1]]/f.w[η_old] *
        f.T0[f.temp_vec[1], η_old]/f.T0[η_old, f.temp_vec[1]], 0.0, f.temp_vec
end

@inline function propose_local(mc::DQMC, model::Model, 
    f::Discrete_MBF2_symm, i::Int, slice::Int)

    rand1=rand(1:15)
    px=mod1(f.conf[1,i, slice] + mod(rand1, 4),4)
    py=mod1(f.conf[2,i, slice] + div(rand1, 4),4)
    f.temp_vec[:]=@inbounds f.choices[px, py]

    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return abs2(detratio) * (f.w[f.temp_vec[1]]*f.w[f.temp_vec[2]])/
        (f.w[f.conf[1,i, slice]]*f.w[f.conf[2,i, slice]]), 
            0.0, f.temp_vec
end


@bm function propose_global_from_conf(mc::DQMC, m::Model, 
    f::Union{Discrete_MBF1_symm, Abstract_DiscreteMBF1_X_symm, Discrete_MBF2_symm})
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
        #If the Update is any shuffle, then the factor fac=1
    end

    return abs2(detratio) * fac, 0.0, nothing
end
###############################################
## interactin matrix exponentials (discrete fields)
##############################################

# exp(-power*δτ*V)
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel, 
    f::Union{Discrete_MBF1_symm, Abstract_DiscreteMBF1_X_symm},
    η_vec::Vector{Int8}, power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    sh3 = sinh(f.α * f.x[η_vec[1]]) # α=sqrt(δτ*U)
    ch3 = cosh(f.α * f.x[η_vec[1]])

    eVop[1,1] = ch3
    eVop[2,2] = ch3          
    eVop[1,2] = power*sh3
    eVop[2,1] = power*sh3
    #Note that index 2 equals index 3 for the Ising-z, and index 4 for the Ising-x case.
    return nothing
end

# exp(-power*δτ*V)
@inline function interaction_matrix_exp_op!(mc::DQMC, model::TwoBandModel, 
    f::Discrete_MBF2_symm,
    η_vec::Vector{Int8}, power::Float64=1., eVop::Matrix{G}=mc.stack.field_cache.eVop1) where G

    sh1 = sinh(f.α * f.x[η_vec[1]])   # α=sqrt(δτ*U)
    ch1 = cosh(f.α * f.x[η_vec[1]])
    sh2 = sinh(f.α * f.x[η_vec[2]])
    ch2 = cosh(f.α * f.x[η_vec[2]])

    eVop[1,1] = ch1*ch2+im*power*sh1*sh2
    eVop[2,2] = ch1*ch2-im*power*sh1*sh2         
    eVop[1,2] = power*(ch2*sh1 - im*ch1*sh2)
    eVop[2,1] = power*(ch2*sh1 + im*ch1*sh2)
    return nothing
end


# exp(-power*δτ*V) full (4N×4N) matrix
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, 
    f::Union{Discrete_MBF1_symm, Abstract_DiscreteMBF1_X_symm}, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N

        sh3 = sinh(f.α * f.x[f.conf[1,i,slice]])    # α=sqrt(δτ*U)
        ch3 = cosh(f.α * f.x[f.conf[1,i,slice]])
        result[i,i]         = ch3
        result[i+N,i+N]   = ch3
        result[i,i+N]      = power*sh3
        result[i+N,i]      = power*sh3
    end
    return nothing
end

# exp(-power*δτ*V) full (4N×4N) matrix
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, 
    f::Discrete_MBF2_symm, 
    result::Union{Matrix{G},SparseMatrixCSC}, slice::Int, power::Float64) where G
    N = length(lattice(mc))
    #result.= 0.0 #I am assuming result=eV always, and eV is (set and stays) zero everywhere except for the positions below.
    @inbounds for i=1:N
        sh1 = sinh(f.α * f.x[f.conf[1,i,slice]])  # α=sqrt(δτ*U)
        ch1 = cosh(f.α * f.x[f.conf[1,i,slice]])
        sh2 = sinh(f.α * f.x[f.conf[2,i,slice]])
        ch2 = cosh(f.α * f.x[f.conf[2,i,slice]])

        result[i,i]         = ch1*ch2 + im*power*sh1*sh2
        result[i+N,i+N]   = ch1*ch2 - im*power*sh1*sh2
        result[i,i+N]      = power*(ch2*sh1 - im*ch1*sh2)
        result[i+N,i]      = power*(ch2*sh1 + im*ch1*sh2)
    end
    return nothing
end
