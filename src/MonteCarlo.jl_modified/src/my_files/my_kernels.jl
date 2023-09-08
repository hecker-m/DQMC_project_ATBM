###########################
### occupation + energy + charge-density for symmetry-optimized Discrete_MBF1
###########################
@inline Base.@propagate_inbounds function occupation_kernel(i, flv, N, G::_GM{<: Matrix}, 
    ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    shift = N * (flv - 1)
    return 2 - 2real(G.val[i+shift, i+shift])
end
@inline Base.@propagate_inbounds function kinetic_energy_kernel(mc, model, ::Nothing, 
    G::_GM{<: Matrix}, flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    # <T> = \sum Tji * (Iij - Gij) = - \sum Tji * (Gij - Iij)
    T = mc.stack.hopping_matrix
    N = length(lattice(mc))
    output = zero(eltype(G.val))
    @inbounds @fastmath for i in 1:N, j in 1:N
        # using T Hermitian for better cache friendliness (conj(Tᵢⱼ) = Tⱼᵢ)
        output += 2real(conj(T[i, j]) * (I[i, j] - G.val[i, j]))+
                2real(conj(T[i+N, j+N]) * (I[i, j] - G.val[i+N, j+N]))
    end
    output
end


@inline Base.@propagate_inbounds function full_cdc_kernel(
    mc, ::Model, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, 
    flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    i, j = ij
    α, γ = flv
    G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(mc))
    id = I[i, j] * I[G0l.k, G0l.l] * I[α, γ]

    return 4*real(1 - Gll.val[i+(α-1)*N, i+(α-1)*N]) * real(1 - G00.val[j+(γ-1)*N, j+(γ-1)*N]) +
            2*real(Gl0.val[i+(α-1)*N, j+(γ-1)*N]*(id - G0l.val[j+(γ-1)*N, i+(α-1)*N]))
end

###########################
### interaction energy
###########################
"""
    interaction_energy_kernel(mc, ::HubbardModel, ::Nothing, greens_matrices, ::Nothing)

Computes the interaction energy - U ⟨(n↑ - 0.5)(n↓ - 0.5)⟩. Note that this is a 
model specific method.
"""

@inline Base.@propagate_inbounds function interaction_energy_kernel(mc, model::TwoBandModel, 
        ::Nothing, G::_GM{<: Matrix}, flv)
    return interaction_energy_kernel(mc, model, nothing, G, flv, field(mc)) 
end
@inline Base.@propagate_inbounds function interaction_energy_kernel(mc, model::TwoBandModel, ::Nothing, 
    G::_GM{<: Matrix}, flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    N = length(lattice(model))
    output = zero(eltype(G.val))
    outtmp1 = zero(eltype(G.val))
    outtmp2 = zero(eltype(G.val))
    @inbounds @fastmath for i in 1:N
        outtmp1=(2*real(G.val[i, i+N]+G.val[i+N, i]))^2;   
        outtmp2=2*real(G.val[i, i]+G.val[i+N, i+N])-
            4*real(G.val[i, i]*G.val[i+N, i+N])-
            2*(real(G.val[i, i+N])^2-imag(G.val[i, i+N])^2)-
            2*(real(G.val[i+N, i])^2-imag(G.val[i+N, i])^2);
        output-=outtmp1+outtmp2
    end
    return output * model.U
end
@inline Base.@propagate_inbounds function interaction_energy_kernel(mc, model::TwoBandModel, ::Nothing, 
    G::_GM{<: Matrix}, flv, ::Union{Cont_MBF1, Discrete_MBF1})
    N = length(lattice(model))
    output = zero(eltype(G.val))
    outtmp1 = zero(eltype(G.val))
    outtmp2 = zero(eltype(G.val))
    @inbounds @fastmath for i in 1:N
        outtmp1=(-G.val[i, i+2N]-G.val[i+2N, i]+
            G.val[i+N, i+3N]+G.val[i+3N, i+N])^2;   
        outtmp2=G.val[i, i]*(1-2*G.val[i+2N, i+2N])+G.val[i+2N, i+2N]+
                G.val[i+3N, i+3N]*(1-2*G.val[i+N, i+N])+G.val[i+N, i+N]-
                G.val[i, i+2N]^2-G.val[i+2N, i]^2-
                G.val[i+N, i+3N]^2-G.val[i+3N, i+N]^2+
                2*G.val[i, i+3N]*G.val[i+N, i+2N]+
                2*G.val[i+3N, i]*G.val[i+2N, i+N]+
                2*G.val[i+N, i]*G.val[i+2N, i+3N]+
                2*G.val[i, i+N]*G.val[i+3N, i+2N];
        output-=outtmp1+outtmp2
    end
    return output * model.U
end
@inline Base.@propagate_inbounds function interaction_energy_kernel(mc, model::TwoBandModel, ::Nothing, 
    G::_GM{<: Matrix}, flv, ::Discrete_MBF1_X)
    N = length(lattice(model))
    output = zero(eltype(G.val))
    outtmp1 = zero(eltype(G.val))
    outtmp2 = zero(eltype(G.val))
    @inbounds @fastmath for i in 1:N
        outtmp1=(-G.val[i, i + 3N] - G.val[i + N, i + 2N] - 
                G.val[i + 2N, i + N] - G.val[i + 3N, i])^2;   
        outtmp2=-G.val[i, i + 3N]^2 + G.val[i + N, i + N] - G.val[i + N, i + 2N]^2 - 2*G.val[i, i + 2N]*G.val[i + N, i + 3N] - G.val[i + 2N, i + N]^2 + G.val[i + 2N, i + 2N] - 
        2*G.val[i + N, i + N]*G.val[i + 2N, i + 2N] - 2*G.val[i, i + N]*G.val[i + 2N, i + 3N] - G.val[i + 3N, i]^2 - 2*G.val[i + 2N, i]*G.val[i + 3N, i + N] - 
        2*G.val[i + N, i]*G.val[i + 3N, i + 2N] + G.val[i, i]*(1 - 2*G.val[i + 3N, i + 3N]) + G.val[i + 3N, i + 3N];
        output-=outtmp1+outtmp2
    end
    return output * model.U
end
###
@inline Base.@propagate_inbounds function interaction_energy_kernel(mc, model::TwoBandModel, ::Nothing, 
    G::_GM{<: Matrix}, flv, ::Union{Cont_MBF2, Discrete_MBF2})
    N = length(lattice(model))
    output = zero(eltype(G.val))
    outtmp1 = zero(eltype(G.val))
    outtmp2 = zero(eltype(G.val))
    @inbounds @fastmath for i in 1:N
        outtmp1=4*(G.val[i, i+3N]+G.val[i+2N, i+N])*
            (G.val[i+N, i+2N]+G.val[i+3N, i]);  
            
        outtmp2=2*G.val[i, i]+2*(1-2*G.val[i,i])*G.val[i+3N, i+3N]+            
            2*G.val[i+N, i+N]+2*(1-2*G.val[i+N,i+N])*G.val[i+2N, i+2N]-
            4*G.val[i, i+2N]*G.val[i+N, i+3N]-
            4*G.val[i+2N, i]*G.val[i+3N, i+N];

        output-=outtmp1+outtmp2
    end
    return output * model.U
end
###
@inline Base.@propagate_inbounds function interaction_energy_kernel(mc, model::TwoBandModel, ::Nothing, 
        G::_GM{<: Matrix}, flv, ::Union{Cont_MBF3, Discrete_MBF3})
    N = length(lattice(model))
    output = zero(eltype(G.val))
    outtmp1 = zero(eltype(G.val))
    outtmp2 = zero(eltype(G.val))
    @inbounds @fastmath for i in 1:N
        outtmp1=4*(G.val[i, i+3N]+G.val[i+2N, i+N])*
            (G.val[i+N, i+2N]+G.val[i+3N, i])+
            (G.val[i, i+2N]+G.val[i+2N, i]
            -G.val[i+N, i+3N]-G.val[i+3N, i+N])^2;  
            
        outtmp2=G.val[i, i]*(1-2*G.val[i+2N, i+2N])+G.val[i+2N, i+2N]+
            G.val[i+3N, i+3N]*(1-2*G.val[i+N, i+N])+G.val[i+N, i+N]-
            G.val[i, i+2N]^2-G.val[i+2N, i]^2-
            G.val[i+N, i+3N]^2-G.val[i+3N, i+N]^2+
            2*G.val[i, i+3N]*G.val[i+N, i+2N]+
            2*G.val[i+3N, i]*G.val[i+2N, i+N]+
            2*G.val[i+N, i]*G.val[i+2N, i+3N]+
            2*G.val[i, i+N]*G.val[i+3N, i+2N]+
            2*G.val[i, i]+2*(1-2*G.val[i,i])*G.val[i+3N, i+3N]+            
            2*G.val[i+N, i+N]+2*(1-2*G.val[i+N,i+N])*G.val[i+2N, i+2N]-
            4*G.val[i, i+2N]*G.val[i+N, i+3N]-
            4*G.val[i+2N, i]*G.val[i+3N, i+N];

        output-=outtmp1+outtmp2
    end
    return output * model.U
end





###########################
### magnetization
###########################

"""
Calculates the ⟨M^x_z(i,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function Mx_z_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, flv) 
    return Mx_z_kernel(mc, model, i, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function Mx_z_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, 
    flv, ::AbstractMagnBosonField)
    N = length(lattice(model))
    return -G.val[i, i+2N]-G.val[i+2N, i]+ 
        G.val[i+N, i+3N]+G.val[i+3N, i+N]
end
@inline Base.@propagate_inbounds function Mx_z_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, 
    flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    N = length(lattice(model))
    return -2*real(G.val[i, i+N]+G.val[i+N, i])
    #Note that for the Ising-x field, the Mx_z_kernel is actually the Mx_x_kernel.
    #Just too lazy to freshly implement.
end

@inline Base.@propagate_inbounds function Mx_z_kernel_old(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, 
    flv,::AbstractMagnBosonField)
    nflav=model.nflav
    N = length(lattice(model))
    output = zero(eltype(G.val))
    @inbounds @fastmath for α=1:nflav, αp=1:nflav
        output+=mc.stack.Mxs[3,α,αp]*(I[αp, α]-G.val[i+(αp-1)*N, i+(α-1)*N])
    end
    return output
end
"""
Calculates the ⟨M^x_x(i,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function Mx_x_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, flv) 
    return Mx_x_kernel(mc, model, i, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function Mx_x_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, 
    flv, ::AbstractMagnBosonField)
    N = length(lattice(model))
    return -G.val[i, i + 3N] - G.val[i + N, i + 2N] - G.val[i + 2N, i + N] - 
    G.val[i + 3N, i]
end
@inline Base.@propagate_inbounds function Mx_x_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, 
    flv, ::Discrete_MBF1_X_symm)
    N = length(lattice(model))
    return -2*real(G.val[i, i+N]+G.val[i+N, i])
   
end

"""
Calculates the ⟨M^x_y(i,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function Mx_y_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, flv) 
    return Mx_y_kernel(mc, model, i, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function Mx_y_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, 
    flv, ::AbstractMagnBosonField)
    N = length(lattice(model))
    return 1i *(-G.val[i, i + 3N] + G.val[i + N, i + 2N] - G.val[i + 2N, i + N] + 
    G.val[i + 3N, i])
end
@inline Base.@propagate_inbounds function Mx_y_kernel(mc, model::TwoBandModel, i, G::_GM{<: Matrix}, 
    flv, ::Discrete_MBF1_X_symm)
    N = length(lattice(model))
    return 2*imag(G.val[i, i+N]-G.val[i+N, i])
   
end
###########################
### magnetic susceptibility
###########################


"""
Calculates the ⟨M^x_z(i,τ) M^x_z(0,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_sdc_Mx_z_kernel(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_sdc_Mx_z_kernel(mc, model, ij, (G, G, G, G), flv, field)
end
@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel(
    mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, 
    flv, ::Discrete_MBF1_symm)
i, j = ij   #we treat j as j=0
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(model))
id = I[i, j] * I[G0l.k, G0l.l] 

4*real(G00.val[j, j+N]+G00.val[j+N, j]) *real(Gll.val[i, i+N]+Gll.val[i+N, i])+
2*real((id-G0l.val[j+N, i+N]) * Gl0.val[i, j])+
2*real((id-G0l.val[j, i]) * Gl0.val[i+N, j+N])-
2*real(G0l.val[j, i+N]*Gl0.val[i, j+N])-
2*real(G0l.val[j+N, i]*Gl0.val[i+N, j])
end
@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel(
    mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_X_symm
    )
i, j = ij   #we treat j as j=0
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(model))
id = I[i, j] * I[G0l.k, G0l.l] 

2*real((id-G0l.val[j, i]) * conj(Gl0.val[i+N, j+N]))+
2*real((id-G0l.val[j+N, i+N]) * conj(Gl0.val[i, j]))+
2*real(G0l.val[j, i+N] * conj(Gl0.val[i, j+N]))+
2*real(G0l.val[j+N, i] * conj(Gl0.val[i+N, j]))
end


@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel(
        mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
    )
    i, j = ij   #we treat j as j=0
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[i, j] * I[G0l.k, G0l.l] 
   
    (G00.val[j, j+2N]-G00.val[j+N, j+3N]+ G00.val[j+2N, j]-G00.val[j+3N, j+N])*
    (Gll.val[i, i+2N]-Gll.val[i+N, i+3N]+ Gll.val[i+2N, i]-Gll.val[i+3N, i+N])+
    (id-G0l.val[j+2N, i+2N]) * Gl0.val[i, j]+
    (id-G0l.val[j+3N, i+3N]) * Gl0.val[i+N, j+N]+
    (id-G0l.val[j, i]) * Gl0.val[i+2N, j+2N]+
    (id-G0l.val[j+N, i+N]) * Gl0.val[i+3N, j+3N]+
    G0l.val[j+3N, i+2N] * Gl0.val[i, j+N]-
    G0l.val[j, i+2N] * Gl0.val[i, j+2N]+
    G0l.val[j+N, i+2N] * Gl0.val[i, j+3N]+
    G0l.val[j+2N, i+3N] * Gl0.val[i+N, j]+
    G0l.val[j, i+3N] * Gl0.val[i+N, j+2N]-
    G0l.val[j+N, i+3N] * Gl0.val[i+N, j+3N]-
    G0l.val[j+2N, i] * Gl0.val[i+2N, j]+
    G0l.val[j+3N, i] * Gl0.val[i+2N, j+N]+
    G0l.val[j+N, i] * Gl0.val[i+2N, j+3N]+
    G0l.val[j+2N, i+N] * Gl0.val[i+3N, j]-
    G0l.val[j+3N, i+N] * Gl0.val[i+3N, j+N]+
    G0l.val[j, i+N] * Gl0.val[i+3N, j+2N]
end

@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel_old(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, G::GreensMatrix, flv)
    return full_sdc_Mx_z_kernel_old(mc, model, ij, (G, G, G, G), flv)
end

@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel_old(
        mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv
    )
    i, j = ij   #we treat j as j=0
	G00, G0l, Gl0, Gll = packed_greens
    nflav=model.nflav
    N = length(lattice(model))
    output = zero(eltype(G00.val))
    id = I[i, j] * I[G0l.k, G0l.l] 

    @inbounds @fastmath for α=1:nflav, αp=1:nflav, γ=1:nflav, γp=1:nflav
        output+=mc.stack.Mxs[3,α,αp]*Gll.val[i+(αp-1)*N, i+(α-1)*N]*
                mc.stack.Mxs[3,γ,γp]*G00.val[j+(γp-1)*N, j+(γ-1)*N]+
                mc.stack.Mxs[3,α,αp]*mc.stack.Mxs[3,γ,γp]*
                (id* I[α, γp]-G0l.val[j+(γp-1)*N, i+(α-1)*N])*Gl0.val[i+(αp-1)*N, j+(γ-1)*N]
    end
    return output
end

@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel_alt(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, G::GreensMatrix, flv)
    return full_sdc_Mx_z_kernel_alt(mc, model, ij, (G, G, G, G), flv)
end

@inline Base.@propagate_inbounds function full_sdc_Mx_z_kernel_alt(
        mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv
    )
    i, j = ij   #we treat j as j=0
	G00, G0l, Gl0, Gll = packed_greens
    nflav=model.nflav
    N = length(lattice(model))
    output = zero(eltype(G00.val))

    @inbounds @fastmath for α=1:nflav, αp=1:nflav, γ=1:nflav, γp=1:nflav
        output+=mc.stack.Mxs[3,α,αp]*Gll.val[i+(αp-1)*N, i+(α-1)*N]*
                mc.stack.Mxs[3,γ,γp]*G00.val[j+(γp-1)*N, j+(γ-1)*N]+
                mc.stack.Mxs[3,α,αp]*mc.stack.Mxs[3,γ,γp]*
                (I[i,j]*I[γp,α]-G0l.val[j+(γp-1)*N, i+(α-1)*N])*Gl0.val[i+(αp-1)*N, j+(γ-1)*N]
    end
    return output
end

"""
Calculates the ⟨M^x_x(i,τ) M^x_x(0,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function full_sdc_Mx_x_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_sdc_Mx_x_kernel(mc, model, ij, G, flv, field(mc))
end

@inline Base.@propagate_inbounds function full_sdc_Mx_x_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_sdc_Mx_x_kernel(mc, model, ij, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function full_sdc_Mx_x_kernel(
    mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_symm
    )
i, j = ij   #we treat j as j=0
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(model))
id = I[i, j] * I[G0l.k, G0l.l] 

2*real((id-G0l.val[j, i]) * conj(Gl0.val[i+N, j+N]))+
2*real((id-G0l.val[j+N, i+N]) * conj(Gl0.val[i, j]))+
2*real(G0l.val[j, i+N] * conj(Gl0.val[i, j+N]))+
2*real(G0l.val[j+N, i] * conj(Gl0.val[i+N, j]))
end
@inline Base.@propagate_inbounds function full_sdc_Mx_x_kernel(
    mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, 
    flv, ::Discrete_MBF1_X_symm)
i, j = ij   #we treat j as j=0
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(model))
id = I[i, j] * I[G0l.k, G0l.l] 

4*real(G00.val[j, j+N]+G00.val[j+N, j]) *real(Gll.val[i, i+N]+Gll.val[i+N, i])+
2*real((id-G0l.val[j+N, i+N]) * Gl0.val[i, j])+
2*real((id-G0l.val[j, i]) * Gl0.val[i+N, j+N])-
2*real(G0l.val[j, i+N]*Gl0.val[i, j+N])-
2*real(G0l.val[j+N, i]*Gl0.val[i+N, j])
end

@inline Base.@propagate_inbounds function full_sdc_Mx_x_kernel(
        mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
    )
    i, j = ij   #we treat j as j=0
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[i, j] * I[G0l.k, G0l.l] 
   
    (G00.val[j, j+3N]+G00.val[j+N, j+2N]+ G00.val[j+2N, j+N]+G00.val[j+3N, j])*
    (Gll.val[i, i+3N]+Gll.val[i+N, i+2N]+ Gll.val[i+2N, i+N]+Gll.val[i+3N, i])+
    (id-G0l.val[j+3N, i+3N]) * Gl0.val[i, j]+
    (id-G0l.val[j+2N, i+2N]) * Gl0.val[i+N, j+N]+
    (id-G0l.val[j+N, i+N]) * Gl0.val[i+2N, j+2N]+
    (id-G0l.val[j, i]) * Gl0.val[i+3N, j+3N]-
    G0l.val[j+2N, i+3N] * Gl0.val[i, j+N]-
    G0l.val[j+N, i+3N] * Gl0.val[i, j+2N]-
    G0l.val[j, i+3N] * Gl0.val[i, j+3N]-
    G0l.val[j+3N, i+2N] * Gl0.val[i+N, j]-
    G0l.val[j+N, i+2N] * Gl0.val[i+N, j+2N]-
    G0l.val[j, i+2N] * Gl0.val[i+N, j+3N]-
    G0l.val[j+3N, i+N] * Gl0.val[i+2N, j]-
    G0l.val[j+2N, i+N] * Gl0.val[i+2N, j+N]-
    G0l.val[j, i+N] * Gl0.val[i+2N, j+3N]-
    G0l.val[j+3N, i] * Gl0.val[i+3N, j]-
    G0l.val[j+2N, i] * Gl0.val[i+3N, j+N]-
    G0l.val[j+N, i] * Gl0.val[i+3N, j+2N]
end
"""
Calculates the ⟨M^x_y(i,τ) M^x_y(0,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function full_sdc_Mx_y_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_sdc_Mx_y_kernel(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function full_sdc_Mx_y_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return full_sdc_Mx_y_kernel(mc, model, ij, (G, G, G, G), flv, field)
end
@inline Base.@propagate_inbounds function full_sdc_Mx_y_kernel(
    mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField
)
i, j = ij   #we treat j as j=0
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(model))
id = I[i, j] * I[G0l.k, G0l.l] 

-(G00.val[j, j+3N]-G00.val[j+N, j+2N]+ G00.val[j+2N, j+N]-G00.val[j+3N, j])*
(Gll.val[i, i+3N]-Gll.val[i+N, i+2N]+ Gll.val[i+2N, i+N]-Gll.val[i+3N, i])+
(id-G0l.val[j+3N, i+3N]) * Gl0.val[i, j]+
(id-G0l.val[j+2N, i+2N]) * Gl0.val[i+N, j+N]+
(id-G0l.val[j+N, i+N]) * Gl0.val[i+2N, j+2N]+
(id-G0l.val[j, i]) * Gl0.val[i+3N, j+3N]+
G0l.val[j+2N, i+3N] * Gl0.val[i, j+N]-
G0l.val[j+N, i+3N] * Gl0.val[i, j+2N]+
G0l.val[j, i+3N] * Gl0.val[i, j+3N]+
G0l.val[j+3N, i+2N] * Gl0.val[i+N, j]+
G0l.val[j+N, i+2N] * Gl0.val[i+N, j+2N]-
G0l.val[j, i+2N] * Gl0.val[i+N, j+3N]-
G0l.val[j+3N, i+N] * Gl0.val[i+2N, j]+
G0l.val[j+2N, i+N] * Gl0.val[i+2N, j+N]+
G0l.val[j, i+N] * Gl0.val[i+2N, j+3N]+
G0l.val[j+3N, i] * Gl0.val[i+3N, j]-
G0l.val[j+2N, i] * Gl0.val[i+3N, j+N]+
G0l.val[j+N, i] * Gl0.val[i+3N, j+2N]
end

@inline Base.@propagate_inbounds function full_sdc_Mx_y_kernel(
        mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_symm
    )
    i, j = ij   #we treat j as j=0
	G00, G0l, Gl0, Gll = packed_greens
    N = length(lattice(model))
    id = I[i, j] * I[G0l.k, G0l.l] 
   
    2*real((id-G0l.val[j, i]) * conj(Gl0.val[i+N, j+N]))+
    2*real((id-G0l.val[j+N, i+N]) * conj(Gl0.val[i, j]))+
    2*real(G0l.val[j, i+N] * conj(Gl0.val[i, j+N]))+
    2*real(G0l.val[j+N, i] * conj(Gl0.val[i+N, j]))
end
@inline Base.@propagate_inbounds function full_sdc_Mx_y_kernel(
    mc, model::TwoBandModel, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, 
    flv, ::Discrete_MBF1_X_symm)
i, j = ij   #we treat j as j=0
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(model))
id = I[i, j] * I[G0l.k, G0l.l] 

4*imag(G00.val[j, j+N]-G00.val[j+N, j]) *imag(Gll.val[i, i+N]-Gll.val[i+N, i])+
2*real((id-G0l.val[j+N, i+N]) * Gl0.val[i, j])+
2*real((id-G0l.val[j, i]) * Gl0.val[i+N, j+N])+
2*real(G0l.val[j, i+N]*Gl0.val[i, j+N])+
2*real(G0l.val[j+N, i]*Gl0.val[i+N, j])
end

################################################################################
### Full CDC_XX kernel
################################################################################
"""
    full_cdc_XX_kernel(mc, model, site_indices, greens_matrices, flavor_indices)

Computes `⟨nᵢ(τ) nⱼ(0)⟩` for the given indices.
"""
@inline Base.@propagate_inbounds function full_cdc_XX_kernel(mc::DQMC, model, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return full_cdc_XX_kernel(mc, model, ij, G, flv, field(mc))
end


@inline Base.@propagate_inbounds function full_cdc_XX_kernel(mc, model, ij, G::GreensMatrix, flv, field::AbstractField)
    return full_cdc_XX_kernel(mc, model, ij, (G, G, G, G), flv, field)
end

@inline Base.@propagate_inbounds function full_cdc_XX_kernel(
    mc, ::Model, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::AbstractField
)
i, j = ij
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(mc))
id = I[G0l.k, G0l.l]

return -(G0l.val[j + 2N, i + 2N]*Gl0.val[i, j]) - G0l.val[j + 3N, i + 2N]*Gl0.val[i, j + N] - G0l.val[j, i + 2N]*Gl0.val[i, j + 2N] - G0l.val[j + N, i + 2N]*Gl0.val[i, j + 3N] - 
    G0l.val[j + 2N, i + 3N]*Gl0.val[i + N, j] - G0l.val[j + 3N, i + 3N]*Gl0.val[i + N, j + N] - G0l.val[j, i + 3N]*Gl0.val[i + N, j + 2N] - G0l.val[j + N, i + 3N]*Gl0.val[i + N, j + 3N] - 
    G0l.val[j + 2N, i]*Gl0.val[i + 2N, j] - G0l.val[j + 3N, i]*Gl0.val[i + 2N, j + N] - G0l.val[j, i]*Gl0.val[i + 2N, j + 2N] - G0l.val[j + N, i]*Gl0.val[i + 2N, j + 3N] - 
    G0l.val[j + 2N, i + N]*Gl0.val[i + 3N, j] - G0l.val[j + 3N, i + N]*Gl0.val[i + 3N, j + N] - G0l.val[j, i + N]*Gl0.val[i + 3N, j + 2N] - G0l.val[j + N, i + N]*Gl0.val[i + 3N, j + 3N] + 
    (G00.val[j, j + 2N] + G00.val[j + N, j + 3N] + G00.val[j + 2N, j] + G00.val[j + 3N, j + N])*(Gll.val[i, i + 2N] + Gll.val[i + N, i + 3N] + Gll.val[i + 2N, i] + Gll.val[i + 3N, i + N]) + 
    id*(Gl0.val[i, j] + Gl0.val[i + N, j + N] + Gl0.val[i + 2N, j + 2N] + Gl0.val[i + 3N, j + 3N])*I[j, i]      
end
@inline Base.@propagate_inbounds function full_cdc_XX_kernel(
    mc, ::Model, ij::NTuple{2}, packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_X_symm
)
i, j = ij
G00, G0l, Gl0, Gll = packed_greens
N = length(lattice(mc))
id = I[G0l.k, G0l.l]

return  -2*real(conj(Gl0.val[i + N, j + N])*G0l.val[j, i] + conj(Gl0.val[i, j + N])*G0l.val[j, i + N] + conj(Gl0.val[i + N, j])*G0l.val[j + N, i] + 
    conj(Gl0.val[i, j])*G0l.val[j + N, i + N]) + 2*id*I[j, i]*real(Gl0.val[i, j] + Gl0.val[i + N, j + N])   
end

###########################
### (++)-pairing susceptibility
###########################

"""
Calculates the ⟨Δ†_++(i,0) Δ_++(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_swave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_swave_kernel(mc, model, ij, G, flv, field(mc))
end

@inline Base.@propagate_inbounds function pc_swave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_swave_kernel(mc, model, ij, (G, G, G, G), flv, field)
end

"""
Calculates the ⟨Δ†_++(i,τ) Δ_++(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_swave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField)
    i, j = ij   
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[i, j] * I[G0l.k, G0l.l]
    
    return -4*(-(id-G0l.val[j, i])*(id-G0l.val[j+N, i+N])-
        (id-G0l.val[j+2N, i+2N])*(id-G0l.val[j+3N, i+3N])+
        G0l.val[j, i+N]*G0l.val[j+N, i]+
        G0l.val[j, i+3N]*G0l.val[j+N, i+2N]-
        G0l.val[j, i+2N]*G0l.val[j+N, i+3N]+
        G0l.val[j+2N, i+N]*G0l.val[j+3N, i]-
        G0l.val[j+2N, i]*G0l.val[j+3N, i+N]+
        G0l.val[j+2N, i+3N]*G0l.val[j+3N, i+2N])  
        #the overall minus takes the missing im^2 into account from syt
end
@inline Base.@propagate_inbounds function pc_swave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    i, j = ij   
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[i, j] * I[G0l.k, G0l.l]
    
    return -4*(2id*(-id+real(G0l.val[j, i]+G0l.val[j+N, i+N]))+
        abs2(G0l.val[j, i+N])+abs2(G0l.val[j+N, i])-
        abs2(G0l.val[j, i])- abs2(G0l.val[j+N, i+N]))
        #the overall minus takes the missing im^2 into account from syt
end
"""
Calculates the ⟨Δ_++(i,τ) Δ†_++(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_swave_kernel_conj(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_swave_kernel_conj(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function pc_swave_kernel_conj(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_swave_kernel_conj(mc, model, ij, (G, G, G, G), flv, field)
end
@inline Base.@propagate_inbounds function pc_swave_kernel_conj(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField)
    i, j = ij   
    Gl0= packed_greens[3]
    N = length(lattice(model))   
    return -4*(-Gl0.val[i,j]*Gl0.val[i+N,j+N]-
        Gl0.val[i+2N,j+2N]*Gl0.val[i+3N,j+3N]+
        Gl0.val[i+N,j]*Gl0.val[i,j+N]+
        Gl0.val[i+3N,j]*Gl0.val[i+2N,j+N]-
        Gl0.val[i+2N,j]*Gl0.val[i+3N,j+N]+
        Gl0.val[i+N,j+2N]*Gl0.val[i,j+3N]-
        Gl0.val[i,j+2N]*Gl0.val[i+N,j+3N]+
        Gl0.val[i+3N,j+2N]*Gl0.val[i+2N,j+3N])  
        #the overall minus takes the missing im^2 into account from syt
end
@inline Base.@propagate_inbounds function pc_swave_kernel_conj(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    i, j = ij   
    Gl0= packed_greens[3]
    N = length(lattice(model))   
    return -4*(abs2(Gl0.val[i, j+N])+abs2(Gl0.val[i+N, j])-
        abs2(Gl0.val[i, j])- abs2(Gl0.val[i+N, j+N]))
        #the overall minus takes the missing im^2 into account from syt
end
@inline Base.@propagate_inbounds function pc_swave_kernel_old(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv,::AbstractMagnBosonField)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    nflav=model.nflav
    N = length(lattice(model))
    output = zero(eltype(G0l.val))
    id = I[i, j] * I[G0l.k, G0l.l]

    @inbounds @fastmath for α=1:nflav, αp=1:nflav, γ=1:nflav, γp=1:nflav
        output+=mc.stack.pnys[1,α,αp]*mc.stack.pnys[1,γ,γp]*
        ((id*I[γp,α]-G0l.val[j+(γp-1)*N, i+(α-1)*N])*
        (id*I[γ,αp]-G0l.val[j+(γ-1)*N, i+(αp-1)*N])-
        (id*I[γp,αp]-G0l.val[j+(γp-1)*N, i+(αp-1)*N])*
        (id*I[γ,α]-G0l.val[j+(γ-1)*N, i+(α-1)*N]))    
        end
    return -output  #minus takes the missing im^2 into account from syt
end

"""
Calculates the 0.5*(⟨Δ†_++(i,τ) Δ_++(j,0)⟩+⟨Δ_++(i,τ) Δ†_++(j,0)⟩) kernel
"""
@inline Base.@propagate_inbounds function pc_swave_kernel_symm(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_swave_kernel_symm(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function pc_swave_kernel_symm(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv,  field::AbstractMagnBosonField)
    return pc_swave_kernel_symm(mc, model, ij, (G, G, G, G), flv, field)
end
@inline Base.@propagate_inbounds function pc_swave_kernel_symm(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, field::AbstractMagnBosonField)
    return 0.5*(pc_swave_kernel(mc, model, ij, packed_greens, flv, field)+
        pc_swave_kernel_conj(mc, model, ij, packed_greens, flv, field))
end


###########################
### (+-)-pairing susceptibility
###########################

"""
Calculates the ⟨Δ†_±(i,0) Δ_±(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_spm_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_spm_wave_kernel(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function pc_spm_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_spm_wave_kernel(mc, model, ij, (G, G, G, G), flv, field)
end

"""
Calculates the ⟨Δ†_±(i,τ) Δ_±(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_spm_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField)
    i, j = ij   #we treat j as j=0, i=target, j=source
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[i, j] * I[G0l.k, G0l.l]

    return -4*(-(id-G0l.val[j, i])*(id-G0l.val[j+N, i+N])-
        (id-G0l.val[j+2N, i+2N])*(id-G0l.val[j+3N, i+3N])+
        G0l.val[j, i+N]*G0l.val[j+N, i]-
        G0l.val[j, i+3N]*G0l.val[j+N, i+2N]+
        G0l.val[j, i+2N]*G0l.val[j+N, i+3N]-
        G0l.val[j+2N, i+N]*G0l.val[j+3N, i]+
        G0l.val[j+2N, i]*G0l.val[j+3N, i+N]+
        G0l.val[j+2N, i+3N]*G0l.val[j+3N, i+2N])  
        #the overall minus takes the missing im^2 into account from syt
end
@inline Base.@propagate_inbounds function pc_spm_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[i, j] * I[G0l.k, G0l.l]
    return -4*(2id*(-id+real(G0l.val[j, i] + G0l.val[j+N, i+N]))-
        abs2(G0l.val[j, i+N])- abs2(G0l.val[j+N, i])-
        abs2(G0l.val[j, i])- abs2(G0l.val[j+N, i+N]))
        #the overall minus takes the missing im^2 into account from syt
end
"""
Calculates the ⟨Δ_±(i,τ) Δ†_±(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_spm_wave_kernel_conj(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_spm_wave_kernel_conj(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function pc_spm_wave_kernel_conj(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_spm_wave_kernel_conj(mc, model, ij, (G, G, G, G), flv, field)
end
@inline Base.@propagate_inbounds function pc_spm_wave_kernel_conj(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField)
    i, j = ij   #we treat j as j=0
    Gl0= packed_greens[3]
    N = length(lattice(model))

    return -4*(-Gl0.val[i,j]*Gl0.val[i+N,j+N]-
        Gl0.val[i+2N,j+2N]*Gl0.val[i+3N,j+3N]+
        Gl0.val[i+N, j]*Gl0.val[i, j+N]-
        Gl0.val[i+3N, j]*Gl0.val[i+2N, j+N]+
        Gl0.val[i+2N, j]*Gl0.val[i+3N, j+N]-
        Gl0.val[i+N, j+2N]*Gl0.val[i, j+3N]+
        Gl0.val[i, j+2N]*Gl0.val[i+N, j+3N]+
        Gl0.val[i+3N, j+2N]*Gl0.val[i+2N, j+3N])  
        #the overall minus takes the missing im^2 into account from syt
end
@inline Base.@propagate_inbounds function pc_spm_wave_kernel_conj(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::Union{Discrete_MBF1_symm, Discrete_MBF1_X_symm})
    i, j = ij   #we treat j as j=0
    Gl0= packed_greens[3]
    N = length(lattice(model))
    return -4*(-abs2(Gl0.val[i, j+N])-abs2(Gl0.val[i+N, j])-
        abs2(Gl0.val[i, j])- abs2(Gl0.val[i+N, j+N]))
        #the overall minus takes the missing im^2 into account from syt
end

@inline Base.@propagate_inbounds function pc_spm_wave_kernel_old(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv,::AbstractMagnBosonField)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    nflav=model.nflav
    N = length(lattice(model))
    output = zero(eltype(G0l.val))
    id = I[i, j] * I[G0l.k, G0l.l]

    @inbounds @fastmath for α=1:nflav, αp=1:nflav, γ=1:nflav, γp=1:nflav
        output+=mc.stack.pnys[2,α,αp]*mc.stack.pnys[2,γ,γp]*
        ((id*I[γp,α]-G0l.val[j+(γp-1)*N, i+(α-1)*N])*(id*I[γ,αp]-G0l.val[j+(γ-1)*N, i+(αp-1)*N])-
        (id*I[γp,αp]-G0l.val[j+(γp-1)*N, i+(αp-1)*N])*(id*I[γ,α]-G0l.val[j+(γ-1)*N, i+(α-1)*N]))    
        end
    return -output  #minus takes the missing im^2 into account from syt
end


"""
Calculates the 0.5*(⟨Δ†_±(i,τ) Δ_±(j,0)⟩+⟨Δ_±(i,τ) Δ†_±(j,0)⟩) kernel
"""
@inline Base.@propagate_inbounds function pc_spm_wave_kernel_symm(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_spm_wave_kernel_symm(mc, model, ij, G, flv, field(mc))
end

@inline Base.@propagate_inbounds function pc_spm_wave_kernel_symm(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_spm_wave_kernel_symm(mc, model, ij, (G, G, G, G), flv, field)
end
@inline Base.@propagate_inbounds function pc_spm_wave_kernel_symm(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, field::AbstractMagnBosonField)
    return 0.5*(pc_spm_wave_kernel(mc, model, ij, packed_greens, flv, field)+
        pc_spm_wave_kernel_conj(mc, model, ij, packed_greens, flv, field))
end

###########################
### XX-pairing susceptibility
###########################

"""
Calculates the ⟨Δ†_X(i,0) Δ_X(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_XX_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_XX_wave_kernel(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function pc_XX_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_XX_wave_kernel(mc, model, ij, (G, G, G, G), flv, field)
end

"""
Calculates the ⟨Δ†_X(i,τ) Δ_X(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_XX_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[G0l.k, G0l.l]

    return 4*(G0l.val[j + N, i + 3N]*G0l.val[j + 2N, i] - G0l.val[j + N, i + 2N]*G0l.val[j + 2N, i + N] - G0l.val[j + N, i]*G0l.val[j + 2N, i + 3N] - 
    G0l.val[j, i + 3N]*G0l.val[j + 3N, i] + G0l.val[j, i + 2N]*G0l.val[j + 3N, i + N] - G0l.val[j, i + N]*G0l.val[j + 3N, i + 2N] + 
    (G0l.val[j + N, i + N] - id*I[j, i])*(G0l.val[j + 2N, i + 2N] - id*I[j, i]) + (G0l.val[j, i] - id*I[j, i])*(G0l.val[j + 3N, i + 3N] - id*I[j, i])) 

end
@inline Base.@propagate_inbounds function pc_XX_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_X_symm)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[G0l.k, G0l.l]
    return 8*(id*I[j, i]*(id*I[j, i] - real(G0l.val[j, i] + G0l.val[j + N, i + N])) + 
        real(-(G0l.val[j, i + N]*G0l.val[j + N, i]) + G0l.val[j, i]*G0l.val[j + N, i + N]))
end

###########################
### YY-onsite triplet pairing susceptibility
###########################
"""
Calculates the ⟨Δ†_Y(i,0) Δ_Y(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_YYzz_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_YYzz_wave_kernel(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function pc_YYzz_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_YYzz_wave_kernel(mc, model, ij, (G, G, G, G), flv, field)
end

"""
Calculates the ⟨Δ†_Y(i,τ) Δ_Y(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_YYzz_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[G0l.k, G0l.l]

    return 4*(-(G0l.val[j, i + 2N]*G0l.val[j + 2N, i]) + G0l.val[j, i + 3N]*G0l.val[j + 2N, i + N] - G0l.val[j, i + N]*G0l.val[j + 2N, i + 3N] + G0l.val[j + N, i + 2N]*G0l.val[j + 3N, i] - 
    G0l.val[j + N, i + 3N]*G0l.val[j + 3N, i + N] - G0l.val[j + N, i]*G0l.val[j + 3N, i + 2N] + (G0l.val[j, i] - id*I[j, i])*(G0l.val[j + 2N, i + 2N] - id*I[j, i]) + 
    (G0l.val[j + N, i + N] - id*I[j, i])*(G0l.val[j + 3N, i + 3N] - id*I[j, i]))
end
@inline Base.@propagate_inbounds function pc_YYzz_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_X_symm)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[G0l.k, G0l.l]
    return 8*(id*I[j, i]*(id*I[j, i] - real(G0l.val[j, i] + G0l.val[j + N, i + N])) + 
        real(conj(G0l.val[j + N, i])*G0l.val[j, i + N] + conj(G0l.val[j, i])*G0l.val[j + N, i + N]))
end

"""
Calculates the ∑_{z,0,x}  ⟨Δ†_Y(i,0) Δ_Y(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_YYsum_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::Union{GreensMatrix, _GM4{T}}, flv) where {T <: Matrix}
    return pc_YYsum_wave_kernel(mc, model, ij, G, flv, field(mc))
end
@inline Base.@propagate_inbounds function pc_YYsum_wave_kernel(mc::DQMC, model::TwoBandModel, ij::NTuple{2}, 
    G::GreensMatrix, flv, field::AbstractMagnBosonField)
    return pc_YYsum_wave_kernel(mc, model, ij, (G, G, G, G), flv, field)
end

"""
Calculates the ∑_{z,0,x} ⟨Δ†_Y(i,τ) Δ_Y(j,0)⟩ kernel
"""
@inline Base.@propagate_inbounds function pc_YYsum_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::AbstractMagnBosonField)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[G0l.k, G0l.l]

    return -4*(G0l.val[j + N, i + 2N]*G0l.val[j + 2N, i + N] - G0l.val[j + N, i]*G0l.val[j + 2N, i + 3N] + G0l.val[j, i + 3N]*G0l.val[j + 3N, i] + 
    G0l.val[j, i + 2N]*(2*G0l.val[j + 2N, i] + G0l.val[j + 3N, i + N]) + G0l.val[j + N, i + 3N]*(G0l.val[j + 2N, i] + 2*G0l.val[j + 3N, i + N]) - 
    G0l.val[j, i + N]*G0l.val[j + 3N, i + 2N] - G0l.val[j, i]*(2*G0l.val[j + 2N, i + 2N] + G0l.val[j + 3N, i + 3N] - 3*id*I[j, i]) - 
    G0l.val[j + N, i + N]*(G0l.val[j + 2N, i + 2N] + 2*G0l.val[j + 3N, i + 3N] - 3*id*I[j, i]) + 
    3*id*I[j, i]*(G0l.val[j + 2N, i + 2N] + G0l.val[j + 3N, i + 3N] - 2*id*I[j, i]))
end
@inline Base.@propagate_inbounds function pc_YYsum_wave_kernel(mc, model, ij::NTuple{2}, 
    packed_greens::_GM4{<: Matrix}, flv, ::Discrete_MBF1_X_symm)
    i, j = ij   #we treat j as j=0
    G0l= packed_greens[2]
    N = length(lattice(model))
    id = I[G0l.k, G0l.l]
    return 8*(3*id*I[j, i]*(id*I[j, i] - real(G0l.val[j, i] + G0l.val[j + N, i + N])) + 
    real(-(G0l.val[j, i + N]*G0l.val[j + N, i]) + (2*conj(G0l.val[j, i]) + G0l.val[j, i])*G0l.val[j + N, i + N]))
end

###########################
### current-current correlations
###########################


@inline Base.@propagate_inbounds function my_cc_kernel(mc, m::Model, sites, G::_GM, flv)
    return my_cc_kernel(mc, m, sites, (G, G, G, G), flv)
end
    
    
    # Basic full Matrix
@inline Base.@propagate_inbounds function my_cc_kernel(
        mc, ::Model, sites::NTuple{4, Int}, 
        packed_greens::_GM4{<: Matrix}, 
        flv::NTuple{2, Int})
    
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
    