


unique_flavors(::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm}) = 2
interaction_eltype(f::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm}) =Float64 

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
function pad_to_unique_flavors(f::Discrete_MBF1_X_symm, m::TwoBandModel, mat)
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
    f::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm}, i::Int, slice::Int)
    η_old = f.conf[1,i, slice]

    f.temp_vec[:]=@inbounds [f.choices[η_old, rand(1:3)]]

    detratio = calculate_detratio!(mc, model , i, f.temp_vec)
    return abs2(detratio) * f.w[f.temp_vec[1]]/f.w[η_old], 0.0, f.temp_vec
end

@bm function propose_global_from_conf(mc::DQMC, m::Model, 
    f::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
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
    f::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm},
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


# exp(-power*δτ*V) full (4N×4N) matrix
@inline function interaction_matrix_exp!(mc::DQMC, model::TwoBandModel, 
    f::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm}, 
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

