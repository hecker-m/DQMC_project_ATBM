"""
Square lattice with next-to-nearest neighbor hopping couplings
"""
function SquareLattice_NNN(Lx, Ly = Lx)
    uc = UnitCell(
        "NNN - Square",
        (Float64[1, 0], Float64[0, 1]),
        [Float64[0, 0]],
        [   
            #on-site bond, due to chem. potential
            Bond(1, 1, (0,  0), 0),

            #NN bonds plus reversed
            Bond(1, 1, ( 1,  0), 1),
            Bond(1, 1, ( 0,  1), 1),
            Bond(1, 1, (-1,  0), 2),
            Bond(1, 1, ( 0, -1), 2),

            #NNN bonds plus reversed
            Bond(1, 1, ( 1,  1), 3),
            Bond(1, 1, ( -1,  1), 3),
            Bond(1, 1, ( -1,  -1), 4),
            Bond(1, 1, ( 1,  -1), 4)
        ]
    )

    Lattice(uc, (Lx, Ly))
end


"""
    HubbardModel(lattice; params...)
    HubbardModel(L, dims; params...)
    HubbardModel(params::Dict)
    HubbardModel(params::NamedTuple)
    HubbardModel(; params...)

Defines an TwoBandModel model on a given (or derived) `lattice`. If a linear system 
size `L` and dimensionality `dims` is given, the`lattice` will be a Cubic 
lattice of fitting size.

Additional parameters (keyword arguments) include:
* `l::AbstractLattice = lattice`: The lattice the model uses. The keyword 
argument takes precedence over the argument `lattice` as well as `dims` and `L`.
* `U::Float64 = 1.0` is the value of the Hubbard Interaction. A positive `U` 
leads to an attractive Hubbard model, a negative to a repulsive.
* `t::Float64 = 1.0` is the hopping strength.
* `mu::Float64` is the chemical potential. In the repulsive model `mu != 0` 
leads to a sign problem.

The model will select an appropriate Hirsch field based on `U`.
"""
mutable struct TwoBandModel{LT <: AbstractLattice} <: Model
    tx::Vector{Float64}
    ty::Vector{Float64}
    tp::Vector{Float64}
    μs::Vector{Float64}
    U::Float64
    l::LT
    peierls::Bool
    Bfield::Float64
    nflav::Int64
end

@inline function TwoBandModel(tx::Vector{<:Real},ty::Vector{<:Real},tp::Vector{<:Real},
        μs::Vector{<:Real}, U::Real, l::AbstractLattice, peierls::Bool,Bfield::Float64,nflav::Int64)
    TwoBandModel(Vector{Float64}(tx),Vector{Float64}(ty), Vector{Float64}(tp),Vector{Float64}(μs),
        Float64(U), l, peierls, Bfield::Float64,nflav::Int64)
end
# Bfield = n0/L^2 has per default n0=1 flux quanta
function TwoBandModel(; 
        dims = 2, L = 2, l = choose_lattice(TwoBandModel, dims, L), 
        U = 1.0, tx = [0.6; 0.2], ty = [0.6; 0.2],tp = [0.8; 1.0], μs = [2.0;-2.0], 
        peierls=false, Bfield=1*1/L^2, nflav=4
    )
    TwoBandModel(tx, ty, tp, μs, U, l, peierls, Bfield, nflav)
end
TwoBandModel(params::Dict{Symbol}) = TwoBandModel(; params...)
TwoBandModel(params::NamedTuple) = TwoBandModel(; params...)
TwoBandModel(lattice::AbstractLattice; kwargs...) = TwoBandModel(l = lattice; kwargs...)
TwoBandModel(L, dims; kwargs...) = TwoBandModel(dims = dims, L = L; kwargs...)


function choose_lattice(::Type{<: TwoBandModel}, dims, L)
    if dims == 1
        return Chain(L)
    elseif dims == 2
        return SquareLattice_NNN(L)
    else
        return CubicLattice(dims, L)
    end
end

# cosmetics
import Base.summary
import Base.show
function Base.summary(model::TwoBandModel)
    " TwoBand-Model"
end
function Base.show(io::IO, model::TwoBandModel)
    println(io, " TwoBand-Model")
    println(io, "\tParameters tx = ($(model.tx[1]),$(model.tx[2])),ty = ($(model.ty[1]),$(model.ty[2])),tp = ($(model.tp[1]),$(model.tp[2]))
     \t µ = ($(model.μs[1]),$(model.μs[2])), U = $(model.U), Peierls $(model.peierls)")
    print(io, "\t$(nameof(typeof(model.l))) with $(length(model.l)) sites")
end
Base.show(io::IO, m::MIME"text/plain", model::TwoBandModel) = print(io, model)

# Convenience
@inline parameters(m::TwoBandModel) = (N = length(m.l), Ls = size(m.l), tx = m.tx, ty = m.ty, tp = m.tp, 
μ = m.μs, U = m.U)

# implement DQMC interface:
@inline lattice(m::TwoBandModel) = m.l
total_flavors(::TwoBandModel) = 4 # 2 spins + 2 bands
unique_flavors(::TwoBandModel) = 2  #changed to 2!!

hopping_eltype(model::TwoBandModel) = model.peierls ? ComplexF64 : typeof(model.tx[1])
hopping_matrix_type(model::TwoBandModel)= Matrix{hopping_eltype(model)}




δ(x::Int,y::Int)= x == y ? 1 : 0 ;  #Kronecker Delta
pbc=true;  # by default, periodic boundary conditions are on
pbcP(i,L)=pbc && i==L+1 ? 1 : i;
pbcM(i,L)=pbc && i==0 ? L : i;

# hopping kernel without Peierls phases
MhopWO(m::TwoBandModel,ix,iy,jx,jy,γ,ν,γp,νp,L)=m.tx[γ]*δ(γ,γp)*δ(ν,νp)*δ(iy,jy)*
    (δ(ix,pbcP(jx+1,L))+δ(ix,pbcM(jx-1,L)))+
    m.ty[γ]*δ(γ,γp)*δ(ν,νp)*δ(ix,jx)*(δ(iy,pbcP(jy+1,L))+δ(iy,pbcM(jy-1,L)))+ 
    m.tp[γ]*δ(γ,γp)*δ(ν,νp)*
    (δ(ix,pbcP(jx+1,L))*δ(iy,pbcP(jy+1,L))+
    δ(ix,pbcM(jx-1,L))*δ(iy,pbcM(jy-1,L))+
    δ(ix,pbcP(jx+1,L))*δ(iy,pbcM(jy-1,L))+
    δ(ix,pbcM(jx-1,L))*δ(iy,pbcP(jy+1,L)));

# hopping kernel with Peierls phases
MhopPeierls(m::TwoBandModel,ix,iy,jx,jy,γ,ν,γp,νp,L,B)=m.tx[γ]*δ(γ,γp)*δ(ν,νp)*δ(iy,jy)*(δ(ix,pbcP(jx+1,L))*exp(im*B*2*pi*jy)+
            δ(ix,pbcM(jx-1,L))*exp(-im*B*2*pi*jy))+
            m.ty[γ]*δ(γ,γp)*δ(ν,νp)*δ(ix,jx)*(δ(iy,pbcP(jy+1,L))*((1-δ(jy,L))+δ(jy,L)*exp(-im*B*2*pi*L*ix))+
            δ(iy,pbcM(jy-1,L))*((1-δ(jy,1))+δ(jy,1)*exp(im*B*2*pi*L*ix)))+
            m.tp[γ]*δ(γ,γp)*δ(ν,νp)*
            (δ(ix,pbcP(jx+1,L))*δ(iy,pbcP(jy+1,L))*((1-δ(jy,L))*exp(im*B*2*pi*(jy+0.5))+
            δ(jy,L)*exp(-im*B*2*pi*(L*jx-0.5)))+
            δ(ix,pbcM(jx-1,L))*δ(iy,pbcM(jy-1,L))*((1-δ(jy,1))*exp(-im*B*2*pi*(jy-0.5))+
            δ(jy,1)*exp(im*B*2*pi*(L*jx-L-0.5)))+
            δ(ix,pbcP(jx+1,L))*δ(iy,pbcM(jy-1,L))*((1-δ(jy,1))*exp(im*B*2*pi*(jy-0.5))+
            δ(jy,1)*exp(im*B*2*pi*(L*jx+L+0.5)))+
            δ(ix,pbcM(jx-1,L))*δ(iy,pbcP(jy+1,L))*((1-δ(jy,L))*exp(-im*B*2*pi*(jy+0.5))+
            δ(jy,L)*exp(-im*B*2*pi*(L*jx+0.5))));
# basis used cᵣ₁↑, cᵣ₁↓, cᵣ₂↑, cᵣ₂↓, with ν denoting spin and γ denoting band
idxM(ix,iy,γ,ν,L)=ix+(iy-1)*L+(ν-1)*(L^2)+(γ-1)*(2*L^2);


"""
    hopping_matrix(model)

Calculates the hopping matrix \$T_{i, j}\$ where \$i, j\$ are
site indices.

This isn't a performance critical method as it is only used once before the
actual simulation.
"""


function hopping_matrix(m::TwoBandModel)
    nflav=m.nflav;
    L=size(m.l)[1];
    N = length(m.l)
    T = zeros(hopping_eltype(m),N*nflav,N*nflav);
    for i=1:2N
        T[i,i]=-m.μs[1];
        T[2N+i,2N+i]=-m.μs[2];
    end
    if m.peierls
        Bvec=m.Bfield*[1; -1; -1; 1];
        for iγ=1:2, iν=1:2, ix = 1:L, iy = 1:L, jγ=1:2, jν=1:2, jx = 1:L, jy = 1:L
            T[idxM(ix,iy,iγ,iν,L),idxM(jx,jy,jγ,jν,L)]+=MhopPeierls(m,ix,iy,jx,jy,iγ,iν,jγ,jν,L,Bvec[iν+2*(iγ-1)]);
        end
    else
        for iγ=1:2, iν=1:2, ix = 1:L, iy = 1:L, jγ=1:2, jν=1:2, jx = 1:L, jy = 1:L
            T[idxM(ix,iy,iγ,iν,L),idxM(jx,jy,jγ,jν,L)]+=MhopWO(m,ix,iy,jx,jy,iγ,iν,jγ,jν,L);
        end
    end
    return m.peierls ? 0.5*(T+T') : T
end



function _save(file::FileLike, entryname::String, m::TwoBandModel)
    write(file, entryname * "/VERSION", 1)
    write(file, entryname * "/tag", "TwoBandModel")
    write(file, entryname * "/tx", m.tx)
    write(file, entryname * "/ty", m.ty)
    write(file, entryname * "/tp", m.tp)
    write(file, entryname * "/mu", m.μs)
    write(file, entryname * "/U", m.U)
    _save(file, entryname * "/l", m.l)
    write(file, entryname * "/peierls", Int(m.peierls))
    write(file, entryname * "/Bfield", m.Bfield)
    write(file, entryname * "/nflav", m.nflav)

    nothing
end

# compat
_load(data, ::Val{:TwoBandModel})=load_model(data, Val(:TwoBandModel));
function load_model(data, ::Val{:TwoBandModel})
    l = _load(data["l"], to_tag(data["l"]))
    TwoBandModel(data["tx"], data["ty"], data["tp"], data["mu"], 
            data["U"], l, Bool(data["peierls"]), data["Bfield"], data["nflav"])
end
#load_model(data, ::Val{:TwoBandModel}) = load_model(data, Val(:TwoBandModel))
field_hint(m, ::Val) = choose_field(m)


