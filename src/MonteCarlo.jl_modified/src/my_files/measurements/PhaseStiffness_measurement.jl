
###########################
### Phase stiffness measurement directly from the current-current susceptibility  
### 
###########################
function phase_stiffness(
    dqmc::DQMC, model::Model; 
        greens_iterator=TimeIntegral(dqmc),
        lattice_iterator = PS_EachBondPairByBravaisDistance(dqmc), wrapper = nothing, 
        flavor_iterator = FlavorIterator(dqmc, 0),#returns constant number of fermion flavors, i.e. 4
        kernel = my_cc_kernel,
        capacity = _default_capacity(dqmc), eltype = geltype(dqmc),
        temp = [zeros(eltype, output_size(lattice_iterator, lattice(dqmc))), 
                zeros(eltype, lattice(dqmc).Ls), zeros(eltype, 4)],
        obs = LogBinner(zeros(eltype, 4), capacity=capacity),
        kwargs...
    )
    li = wrapper === nothing ? lattice_iterator : wrapper(lattice_iterator)
    return DQMCMeasurement(dqmc, model, greens_iterator, li, flavor_iterator, kernel; 
        temp=temp, obs=obs, kwargs...)
end






###########################
### my current-current kernel
###########################

@inline Base.@propagate_inbounds function my_cc_kernel(mc, m::Model, sites, G::Union{GreensMatrix, _GM4{T}}, flv) where T
    return my_cc_kernel(mc, m, sites, G, flv, field(mc))
end

@inline Base.@propagate_inbounds function my_cc_kernel(mc, m::Model, sites, G::_GM, flv, field::AbstractMagnBosonField)
    return my_cc_kernel(mc, m, sites, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function my_cc_kernel(
    mc, ::Model, sites::NTuple{4, Int}, 
    packed_greens::_GM4{<: Matrix}, 
    flv, ::Union{Abstract_DiscreteMBF1_X_symm, Cont_MBF1_X_symm, Discrete_MBF2_symm})

    #(s1, t1, s2, t2)= (r′, r′+ b1, r =r′+d, r+b2)
    #src1, trg1, src2, trg2 = sites      #src1 = r′+b′, trg1 = r′, src2 =r+b, trg2= r
    rp, rpPbp, r, rPb = sites
    G00, G0l, Gl0, Gll = packed_greens
    T = mc.stack.hopping_matrix
    id = I[G0l.k, G0l.l]

    N = length(lattice(mc))


    return -4*real(Gll.val[rPb, r]*T[r, rPb] + Gll.val[N + rPb, N + r]*T[N + r, N + rPb] - Gll.val[r, rPb]*T[rPb, r] - Gll.val[N + r, N + rPb]*T[N + rPb, N + r])*
    real(G00.val[rpPbp, rp]*T[rp, rpPbp] + G00.val[N + rpPbp, N + rp]*T[N + rp, N + rpPbp] - G00.val[rp, rpPbp]*T[rpPbp, rp] - G00.val[N + rp, N + rpPbp]*T[N + rpPbp, N + rp]) + 
   2*real(Gl0.val[rPb, rp]*(G0l.val[rpPbp, r] - id*I[rpPbp, r])*T[r, rPb]*T[rp, rpPbp] + G0l.val[rpPbp, N + r]*Gl0.val[N + rPb, rp]*T[N + r, N + rPb]*T[rp, rpPbp] + 
      G0l.val[N + rpPbp, r]*Gl0.val[rPb, N + rp]*T[r, rPb]*T[N + rp, N + rpPbp] + Gl0.val[N + rPb, N + rp]*(G0l.val[N + rpPbp, N + r] - id*I[rpPbp, r])*T[N + r, N + rPb]*
       T[N + rp, N + rpPbp] + Gl0.val[r, rp]*(-G0l.val[rpPbp, rPb] + id*I[rpPbp, rPb])*T[rp, rpPbp]*T[rPb, r] - G0l.val[N + rpPbp, rPb]*Gl0.val[r, N + rp]*T[N + rp, N + rpPbp]*
       T[rPb, r] - G0l.val[rpPbp, N + rPb]*Gl0.val[N + r, rp]*T[rp, rpPbp]*T[N + rPb, N + r] + Gl0.val[N + r, N + rp]*(-G0l.val[N + rpPbp, N + rPb] + id*I[rpPbp, rPb])*
       T[N + rp, N + rpPbp]*T[N + rPb, N + r] + Gl0.val[rPb, rpPbp]*(-G0l.val[rp, r] + id*I[rp, r])*T[r, rPb]*T[rpPbp, rp] - 
      G0l.val[rp, N + r]*Gl0.val[N + rPb, rpPbp]*T[N + r, N + rPb]*T[rpPbp, rp] + Gl0.val[r, rpPbp]*(G0l.val[rp, rPb] - id*I[rp, rPb])*T[rPb, r]*T[rpPbp, rp] + 
      G0l.val[rp, N + rPb]*Gl0.val[N + r, rpPbp]*T[N + rPb, N + r]*T[rpPbp, rp] - G0l.val[N + rp, r]*Gl0.val[rPb, N + rpPbp]*T[r, rPb]*T[N + rpPbp, N + rp] + 
      Gl0.val[N + rPb, N + rpPbp]*(-G0l.val[N + rp, N + r] + id*I[rp, r])*T[N + r, N + rPb]*T[N + rpPbp, N + rp] + 
      G0l.val[N + rp, rPb]*Gl0.val[r, N + rpPbp]*T[rPb, r]*T[N + rpPbp, N + rp] + Gl0.val[N + r, N + rpPbp]*(G0l.val[N + rp, N + rPb] - id*I[rp, rPb])*T[N + rPb, N + r]*
       T[N + rpPbp, N + rp])
end


@inline Base.@propagate_inbounds function my_cc_kernel(
        mc, ::Model, sites::NTuple{4, Int}, 
        packed_greens::_GM4{<: Matrix}, 
        flv, ::AbstractMagnBosonField)

    #(s1, t1, s2, t2)= (r′, r′+ b1, r =r′+d, r+b2)
    #src1, trg1, src2, trg2 = sites      #src1 = r′+b′, trg1 = r′, src2 =r+b, trg2= r
    rp, rpPbp, r, rPb = sites
    G00, G0l, Gl0, Gll = packed_greens
    T = mc.stack.hopping_matrix
    id = I[G0l.k, G0l.l]

    N = length(lattice(mc))


    return -((Gll.val[rPb, r]*T[r, rPb] + Gll.val[N + rPb, N + r]*T[N + r, N + rPb] + Gll.val[2N + rPb, 2N + r]*T[2N + r, 2N + rPb] + 
        Gll.val[3N + rPb, 3N + r]*T[3N + r, 3N + rPb] - Gll.val[r, rPb]*T[rPb, r] - Gll.val[N + r, N + rPb]*T[N + rPb, N + r] - 
        Gll.val[2N + r, 2N + rPb]*T[2N + rPb, 2N + r] - Gll.val[3N + r, 3N + rPb]*T[3N + rPb, 3N + r])*
        (G00.val[rpPbp, rp]*T[rp, rpPbp] + G00.val[N + rpPbp, N + rp]*T[N + rp, N + rpPbp] + G00.val[2N + rpPbp, 2N + rp]*T[2N + rp, 2N + rpPbp] + 
        G00.val[3N + rpPbp, 3N + rp]*T[3N + rp, 3N + rpPbp] - G00.val[rp, rpPbp]*T[rpPbp, rp] - G00.val[N + rp, N + rpPbp]*T[N + rpPbp, N + rp] - 
        G00.val[2N + rp, 2N + rpPbp]*T[2N + rpPbp, 2N + rp] - G00.val[3N + rp, 3N + rpPbp]*T[3N + rpPbp, 3N + rp]))+
        Gl0.val[rPb, rp]*(G0l.val[rpPbp, r] - id*I[rpPbp, r])*T[r, rPb]*T[rp, rpPbp] + G0l.val[rpPbp, N + r]*Gl0.val[N + rPb, rp]*T[N + r, N + rPb]*T[rp, rpPbp] + 
        G0l.val[rpPbp, 2N + r]*Gl0.val[2N + rPb, rp]*T[2N + r, 2N + rPb]*T[rp, rpPbp] + G0l.val[rpPbp, 3N + r]*Gl0.val[3N + rPb, rp]*T[3N + r, 3N + rPb]*T[rp, rpPbp] + 
        G0l.val[N + rpPbp, r]*Gl0.val[rPb, N + rp]*T[r, rPb]*T[N + rp, N + rpPbp] + Gl0.val[N + rPb, N + rp]*(G0l.val[N + rpPbp, N + r] - id*I[rpPbp, r])*T[N + r, N + rPb]*
        T[N + rp, N + rpPbp] + G0l.val[N + rpPbp, 2N + r]*Gl0.val[2N + rPb, N + rp]*T[2N + r, 2N + rPb]*T[N + rp, N + rpPbp] + 
        G0l.val[N + rpPbp, 3N + r]*Gl0.val[3N + rPb, N + rp]*T[3N + r, 3N + rPb]*T[N + rp, N + rpPbp] + 
        G0l.val[2N + rpPbp, r]*Gl0.val[rPb, 2N + rp]*T[r, rPb]*T[2N + rp, 2N + rpPbp] + G0l.val[2N + rpPbp, N + r]*Gl0.val[N + rPb, 2N + rp]*T[N + r, N + rPb]*
        T[2N + rp, 2N + rpPbp] + Gl0.val[2N + rPb, 2N + rp]*(G0l.val[2N + rpPbp, 2N + r] - id*I[rpPbp, r])*T[2N + r, 2N + rPb]*T[2N + rp, 2N + rpPbp] + 
        G0l.val[2N + rpPbp, 3N + r]*Gl0.val[3N + rPb, 2N + rp]*T[3N + r, 3N + rPb]*T[2N + rp, 2N + rpPbp] + 
        G0l.val[3N + rpPbp, r]*Gl0.val[rPb, 3N + rp]*T[r, rPb]*T[3N + rp, 3N + rpPbp] + G0l.val[3N + rpPbp, N + r]*Gl0.val[N + rPb, 3N + rp]*T[N + r, N + rPb]*
        T[3N + rp, 3N + rpPbp] + G0l.val[3N + rpPbp, 2N + r]*Gl0.val[2N + rPb, 3N + rp]*T[2N + r, 2N + rPb]*T[3N + rp, 3N + rpPbp] + 
        Gl0.val[3N + rPb, 3N + rp]*(G0l.val[3N + rpPbp, 3N + r] - id*I[rpPbp, r])*T[3N + r, 3N + rPb]*T[3N + rp, 3N + rpPbp] + 
        Gl0.val[r, rp]*(-G0l.val[rpPbp, rPb] + id*I[rpPbp, rPb])*T[rp, rpPbp]*T[rPb, r] - G0l.val[N + rpPbp, rPb]*Gl0.val[r, N + rp]*T[N + rp, N + rpPbp]*T[rPb, r] - 
        G0l.val[2N + rpPbp, rPb]*Gl0.val[r, 2N + rp]*T[2N + rp, 2N + rpPbp]*T[rPb, r] - G0l.val[3N + rpPbp, rPb]*Gl0.val[r, 3N + rp]*T[3N + rp, 3N + rpPbp]*T[rPb, r] - 
        G0l.val[rpPbp, N + rPb]*Gl0.val[N + r, rp]*T[rp, rpPbp]*T[N + rPb, N + r] + Gl0.val[N + r, N + rp]*(-G0l.val[N + rpPbp, N + rPb] + id*I[rpPbp, rPb])*T[N + rp, N + rpPbp]*
        T[N + rPb, N + r] - G0l.val[2N + rpPbp, N + rPb]*Gl0.val[N + r, 2N + rp]*T[2N + rp, 2N + rpPbp]*T[N + rPb, N + r] - 
        G0l.val[3N + rpPbp, N + rPb]*Gl0.val[N + r, 3N + rp]*T[3N + rp, 3N + rpPbp]*T[N + rPb, N + r] - 
        G0l.val[rpPbp, 2N + rPb]*Gl0.val[2N + r, rp]*T[rp, rpPbp]*T[2N + rPb, 2N + r] - G0l.val[N + rpPbp, 2N + rPb]*Gl0.val[2N + r, N + rp]*T[N + rp, N + rpPbp]*
        T[2N + rPb, 2N + r] + Gl0.val[2N + r, 2N + rp]*(-G0l.val[2N + rpPbp, 2N + rPb] + id*I[rpPbp, rPb])*T[2N + rp, 2N + rpPbp]*T[2N + rPb, 2N + r] - 
        G0l.val[3N + rpPbp, 2N + rPb]*Gl0.val[2N + r, 3N + rp]*T[3N + rp, 3N + rpPbp]*T[2N + rPb, 2N + r] - 
        G0l.val[rpPbp, 3N + rPb]*Gl0.val[3N + r, rp]*T[rp, rpPbp]*T[3N + rPb, 3N + r] - G0l.val[N + rpPbp, 3N + rPb]*Gl0.val[3N + r, N + rp]*T[N + rp, N + rpPbp]*
        T[3N + rPb, 3N + r] - G0l.val[2N + rpPbp, 3N + rPb]*Gl0.val[3N + r, 2N + rp]*T[2N + rp, 2N + rpPbp]*T[3N + rPb, 3N + r] + 
        Gl0.val[3N + r, 3N + rp]*(-G0l.val[3N + rpPbp, 3N + rPb] + id*I[rpPbp, rPb])*T[3N + rp, 3N + rpPbp]*T[3N + rPb, 3N + r] + 
        Gl0.val[rPb, rpPbp]*(-G0l.val[rp, r] + id*I[rp, r])*T[r, rPb]*T[rpPbp, rp] - G0l.val[rp, N + r]*Gl0.val[N + rPb, rpPbp]*T[N + r, N + rPb]*T[rpPbp, rp] - 
        G0l.val[rp, 2N + r]*Gl0.val[2N + rPb, rpPbp]*T[2N + r, 2N + rPb]*T[rpPbp, rp] - G0l.val[rp, 3N + r]*Gl0.val[3N + rPb, rpPbp]*T[3N + r, 3N + rPb]*T[rpPbp, rp] + 
        Gl0.val[r, rpPbp]*(G0l.val[rp, rPb] - id*I[rp, rPb])*T[rPb, r]*T[rpPbp, rp] + G0l.val[rp, N + rPb]*Gl0.val[N + r, rpPbp]*T[N + rPb, N + r]*T[rpPbp, rp] + 
        G0l.val[rp, 2N + rPb]*Gl0.val[2N + r, rpPbp]*T[2N + rPb, 2N + r]*T[rpPbp, rp] + G0l.val[rp, 3N + rPb]*Gl0.val[3N + r, rpPbp]*T[3N + rPb, 3N + r]*T[rpPbp, rp] - 
        G0l.val[N + rp, r]*Gl0.val[rPb, N + rpPbp]*T[r, rPb]*T[N + rpPbp, N + rp] + Gl0.val[N + rPb, N + rpPbp]*(-G0l.val[N + rp, N + r] + id*I[rp, r])*T[N + r, N + rPb]*
        T[N + rpPbp, N + rp] - G0l.val[N + rp, 2N + r]*Gl0.val[2N + rPb, N + rpPbp]*T[2N + r, 2N + rPb]*T[N + rpPbp, N + rp] - 
        G0l.val[N + rp, 3N + r]*Gl0.val[3N + rPb, N + rpPbp]*T[3N + r, 3N + rPb]*T[N + rpPbp, N + rp] + 
        G0l.val[N + rp, rPb]*Gl0.val[r, N + rpPbp]*T[rPb, r]*T[N + rpPbp, N + rp] + Gl0.val[N + r, N + rpPbp]*(G0l.val[N + rp, N + rPb] - id*I[rp, rPb])*T[N + rPb, N + r]*
        T[N + rpPbp, N + rp] + G0l.val[N + rp, 2N + rPb]*Gl0.val[2N + r, N + rpPbp]*T[2N + rPb, 2N + r]*T[N + rpPbp, N + rp] + 
        G0l.val[N + rp, 3N + rPb]*Gl0.val[3N + r, N + rpPbp]*T[3N + rPb, 3N + r]*T[N + rpPbp, N + rp] - 
        G0l.val[2N + rp, r]*Gl0.val[rPb, 2N + rpPbp]*T[r, rPb]*T[2N + rpPbp, 2N + rp] - G0l.val[2N + rp, N + r]*Gl0.val[N + rPb, 2N + rpPbp]*T[N + r, N + rPb]*
        T[2N + rpPbp, 2N + rp] + Gl0.val[2N + rPb, 2N + rpPbp]*(-G0l.val[2N + rp, 2N + r] + id*I[rp, r])*T[2N + r, 2N + rPb]*T[2N + rpPbp, 2N + rp] - 
        G0l.val[2N + rp, 3N + r]*Gl0.val[3N + rPb, 2N + rpPbp]*T[3N + r, 3N + rPb]*T[2N + rpPbp, 2N + rp] + 
        G0l.val[2N + rp, rPb]*Gl0.val[r, 2N + rpPbp]*T[rPb, r]*T[2N + rpPbp, 2N + rp] + G0l.val[2N + rp, N + rPb]*Gl0.val[N + r, 2N + rpPbp]*T[N + rPb, N + r]*
        T[2N + rpPbp, 2N + rp] + Gl0.val[2N + r, 2N + rpPbp]*(G0l.val[2N + rp, 2N + rPb] - id*I[rp, rPb])*T[2N + rPb, 2N + r]*T[2N + rpPbp, 2N + rp] + 
        G0l.val[2N + rp, 3N + rPb]*Gl0.val[3N + r, 2N + rpPbp]*T[3N + rPb, 3N + r]*T[2N + rpPbp, 2N + rp] + 
        (-(G0l.val[3N + rp, r]*Gl0.val[rPb, 3N + rpPbp]*T[r, rPb]) - G0l.val[3N + rp, N + r]*Gl0.val[N + rPb, 3N + rpPbp]*T[N + r, N + rPb] - 
        G0l.val[3N + rp, 2N + r]*Gl0.val[2N + rPb, 3N + rpPbp]*T[2N + r, 2N + rPb] - G0l.val[3N + rp, 3N + r]*Gl0.val[3N + rPb, 3N + rpPbp]*T[3N + r, 3N + rPb] + 
        id*Gl0.val[3N + rPb, 3N + rpPbp]*I[rp, r]*T[3N + r, 3N + rPb] + G0l.val[3N + rp, rPb]*Gl0.val[r, 3N + rpPbp]*T[rPb, r] + 
        G0l.val[3N + rp, N + rPb]*Gl0.val[N + r, 3N + rpPbp]*T[N + rPb, N + r] + G0l.val[3N + rp, 2N + rPb]*Gl0.val[2N + r, 3N + rpPbp]*T[2N + rPb, 2N + r] + 
        Gl0.val[3N + r, 3N + rpPbp]*(G0l.val[3N + rp, 3N + rPb] - id*I[rp, rPb])*T[3N + rPb, 3N + r])*T[3N + rpPbp, 3N + rp]
end




    # Basic full Matrix
@inline Base.@propagate_inbounds function my_cc_kernel(
        mc, ::Model, sites::NTuple{4, Int}, 
        packed_greens::_GM4{<: Matrix}, 
        flv::NTuple{2, Int}, ::AbstractMagnBosonField)
    
    l, jp, i, j = sites
    G00, G0l, Gl0, Gll = packed_greens
    T = mc.stack.hopping_matrix
    id = I[G0l.k, G0l.l]
    
    N = length(lattice(mc))
    γ, α = flv
    s1 = l + N * (γ - 1); s2 = i + N * (α - 1)
    t1 = jp + N * (γ - 1); t2 = j + N * (α - 1)
    
    T1_st = T[s1, t1]   #Mt[l γ, j' γ]
    T1_ts = conj(T1_st) #Mt[j' γ, l γ]
    T2_st = T[s2, t2]   #Mt[i α, j α]
    T2_ts = conj(T2_st) #Mt[j α, i α]
        #(Mt[l γ, j' γ]*G00.val[j' γ, l γ]-Mt[j' γ, l γ]*  G00.val[l γ, j' γ])
        #*(Mt[j α, i α]*Gll.val[i α, j α] - Mt[i α, j α] * Gll.val[j α, i α])
        #+ (- Mt[j α, i α] * Mt[j' γ, l γ] * (id * I[l γ, j α] - G0l.val[l γ, j α]) * Gl0.val[i α, j' γ] +
        #+ Mt[j α, i α] * Mt[l γ, j' γ] * (id * I[j' γ, j α] - G0l.val[j' γ, j α]) * Gl0.val[i α, l γ] +
        #+ Mt[i α, j α] * Mt[j' γ, l γ] * (id * I[l γ, i α] - G0l.val[l γ, i α]) * Gl0.val[j α, j' γ] +
        #- Mt[i α, j α] * Mt[l γ, j' γ] * (id * I[i α, j' γ] - G0l.val[j' γ, i α]) * Gl0.val[j α, l γ] )
    
        #my notes would require a prefactor 1/4*(j_x-i_x)*(jp_x - l_x)
    output = (
        (T1_st * G00.val[t1, s1] - T1_ts * G00.val[s1, t1])*
        (T2_ts * Gll.val[s2, t2] - T2_st * Gll.val[t2, s2])     
    ) + (
        - T2_ts * T1_ts * (id * I[s1, t2] - G0l.val[s1, t2]) * Gl0.val[s2, t1] +
        + T2_ts * T1_st * (id * I[t1, t2] - G0l.val[t1, t2]) * Gl0.val[s2, s1] +
        + T2_st * T1_ts * (id * I[s1, s2] - G0l.val[s1, s2]) * Gl0.val[t2, t1] +
        - T2_st * T1_st * (id * I[s2, t1] - G0l.val[t1, s2]) * Gl0.val[t2, s1] 
    )
    
    return output
end
    