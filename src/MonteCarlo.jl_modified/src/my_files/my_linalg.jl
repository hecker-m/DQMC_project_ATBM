



@inline function vldiv22!(cache::StandardFieldCache, R::Matrix{G}, Δ::Matrix{G}) where G <:Number
  cache.invRΔ=Δ*inv(R)
  nothing
end
function vsub!(trg::Matrix{ComplexF64}, ::UniformScaling, src::Matrix{ComplexF64}, i::Int, N::Int)
  @inbounds for (α, iα) in enumerate(i:N:size(src, 2)) # α∈{1,..,nflav} and iα∈{2,18,34,50} e.g.
      @inbounds @fastmath for μ in axes(trg, 1) #μ∈{1,...,N*nflav} e.g.          
          trg[μ, α] = - src[μ, iα]
      end
      trg[iα, α] += 1.0
  end
  nothing
end
function vmul!(trg::Matrix{ComplexF64}, M::Matrix{ComplexF64}, src::Matrix{ComplexF64}, i::Int, N::Int)
  r = i:N:size(src, 1)
  @inbounds for μ in axes(trg, 1), α in eachindex(r)
      tij = 0.0
      for (ν, iν) in enumerate(r)  # ν∈{1,..,nflav} and  iν∈{2,18,34,50} e.g.
          tij += M[α, ν] * src[iν,μ]  #M=invRΔ
      end                             
      trg[μ, α] = tij
  end
  nothing
end
# function vsubkron!(G::Matrix{ComplexF64}, L::Matrix{ComplexF64}, R::Matrix{ComplexF64})
#   @inbounds @fastmath for k in axes(L, 1), l in axes(R, 1),  m in axes(L, 2) 
#   #for k in axes(L, 1), l in axes(R, 1), m in axes(L, 2)    
#       G[k, l] -= L[k, m] * R[l, m]
#   end  
# end


@inline vsubkron!(G::AbstractArray, L::AbstractArray,
  R::AbstractArray, LR::AbstractArray) = vsubkron!(G, L, R)

function vsubkron!(G::Matrix{ComplexF64}, L::Matrix{ComplexF64}, R::Matrix{ComplexF64},
    LR::Matrix{ComplexF64})
    vmul!(LR, L, transpose(R))
    #mul!(LR, L, transpose(R))
    G .-= LR
    nothing
end


function udt_AVX_pivot!(
    U::AbstractArray{C, 2}, 
    D::AbstractArray{Float64, 1}, 
    input::AbstractArray{C, 2},
    pivot::AbstractArray{Int64, 1} = Vector(UnitRange(1:size(input, 1))),
    temp::AbstractArray{C, 1} = Vector{C}(undef, length(D)),
    apply_pivot::Val = Val(true)
) where {C <: Complex}

    U[:,:], input[:,:], pivot[:]= @views decompose_udt!(input, D, apply_pivot)
    nothing
end


function decompose_udt!(A::AbstractMatrix{C}, D, ::Val{true}) where C<:Number
    n = length(D)
    F = qr!(A, Val(true))
    R = F.R # F.R is of regular matrix type
  
    @views F.p[F.p] = 1:n
  
    @inbounds for i in 1:n
      D[i] = abs(real(R[i,i]))
    end
  
    lmul!(Diagonal(1 ./ D), R)
  
    return Matrix(F.Q), R[:, F.p], F.p # Q, (D is modified in-place), T  # was full(F.Q) before upgrade
  end

  function decompose_udt!(A::AbstractMatrix{C}, D, ::Val{false}) where C<:Number
    n = length(D)
    F = qr!(A, Val(true))
    R = F.R # F.R is of regular matrix type
  
  
    @inbounds for i in 1:n
      D[i] = abs(real(R[i,i]))
    end
  
    lmul!(Diagonal(1 ./ D), R)
  
    return Matrix(F.Q), R, F.p # Q, (D is modified in-place), T  # was full(F.Q) before upgrade
  end
