
"""
    GlobalConstMove([mc, model])

A global update that shifts entire configuration by a constant.
"""
struct GlobalConstMove <: AbstractGlobalUpdate end
GlobalConstMove(mc, model) = GlobalConstMove()
GlobalConstMove(N) = [GlobalConstMove() for _ in 1:N]
name(::GlobalConstMove) = "GlobalConstMove"
_load(::FileLike, ::Val{:GlobalConstMove}) = GlobalConstMove()


@bm function propose_conf!(::GlobalConstMove, mc, model, field::AbstractContMBF)
    Nϕ=mc.parameters.Nϕ
    c = conf(field); tc = temp_conf(field); δ=field.temp_vec;
    copyto!(tc, c)
    δ[:]=randuniform(mc.parameters.box_global, Nϕ)
    N=length(model.l)
    Nslices=mc.parameters.slices
    # for i in 1:N, ℓ in 1:Nslices
    #     c[:,i,ℓ]+=δ
    # end
    for i in 1:Nϕ
        c[i,:,:] .+= δ[i]
    end
    return nothing
end


"""
    PartialGlobalFlip(p[,mc, model])
A global update of type flip.
A partial global flip that flips a percentage p of the configuration.
"""
mutable struct PartialGlobalFlip <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    PartialGlobalFlip(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
PartialGlobalFlip(p, mc::DQMC, model::Model) = PartialGlobalFlip(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
PartialGlobalFlip(mc::DQMC, model::Model) = PartialGlobalFlip(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
PartialGlobalFlip() = PartialGlobalFlip(1.0,  0, 10000)

name(::PartialGlobalFlip) = "PartialGlobalFlip"
function _save(f::FileLike, name::String, update::PartialGlobalFlip)
    write(f, "$name/tag", :PartialGlobalFlip)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:PartialGlobalFlip}) = PartialGlobalFlip(f["p"], f["count"], f["full_update_interval"])

function init!(mc, u::PartialGlobalFlip)
    resize!(u.indices, length(lattice(mc)))
    u.indices .= 1:length(lattice(mc))
    nothing
end

@bm function propose_conf!(u::PartialGlobalFlip, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    u.count+=1;
    N=length(lattice(mc))
    k=Int8(length(field.x))   #number of discretization (typically k=4)

    if mod(u.count, u.full_update_interval)==0
        Nflipped=N
    else
        Nflipped=max(1, round(Int, u.p *N))
    end
    shuffle!(u.indices)
    for slice in 1:nslices(mc), i in u.indices[1:Nflipped]
        c[:, i, slice] .= Int8(k+1) .- c[:, i, slice]
    end
    return nothing
end

"""
    propose function for GlobalFlip([mc, model])

A global update that flips the configuration (±1 -> ∓1).
"""
@bm function propose_conf!(::GlobalFlip, mc, model, 
    field::Union{Discrete_MBF1, Discrete_MBF1_X, Discrete_MBF1_symm, Discrete_MBF1_X_symm})

    k=Int8(length(field.x))   #number of discretization (typically k=4)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    c .= Int8(k+1) .- c
    return nothing
end

"""
    PartGlobalXorYshift(p[,mc, model])
A global update of type flip.
A partial global shift a number p of rows by tx=a*ex, or a number p of columns by ty=a*ey.
"""
mutable struct pGlobalXorYshift <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    pGlobalXorYshift(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
pGlobalXorYshift(p, mc::DQMC, model::Model) = pGlobalXorYshift(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
pGlobalXorYshift(mc::DQMC, model::Model) = pGlobalXorYshift(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
pGlobalXorYshift() = pGlobalXorYshift(1.0,  0, 10000)

name(::pGlobalXorYshift) = "pGlobalXorYshift"
function _save(f::FileLike, name::String, update::pGlobalXorYshift)
    write(f, "$name/tag", :pGlobalXorYshift)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:pGlobalXorYshift}) = pGlobalXorYshift(f["p"], f["count"], f["full_update_interval"])

function init!(mc, u::pGlobalXorYshift)
    L= mc.model.l.Ls[1]
    resize!(u.indices, L)
    u.indices .= 1:L
    nothing
end

@bm function propose_conf!(u::pGlobalXorYshift, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    L= model.l.Ls[1]
    u.count+=1;
    Bsrcdir2trg = model.l[:Bravais_srcdir2trg]::Matrix{Int}
    N=length(lattice(mc))
    mod(u.count, u.full_update_interval)==0 ? Lshifted=L : Lshifted=max(1, round(Int, u.p *L))

    shuffle!(u.indices)
    if isodd(div(u.count, 3))   #trying three times x, then 3 times y, etc.
            for  row in u.indices[1:Lshifted]
                x_dir=rand([2, L])            #chooses between [+a1 (2) ,-a1 (L)]
                for ix in 1:L
                    src=ix+(row-1)*L
                    trg=Bsrcdir2trg[src, x_dir] # dir=2 is a shift tx=a*ex
                    for slice in 1:nslices(mc)
                        c[:, trg, slice] .= tc[:, src, slice]
                    end
                end
            end
    else
        for  column in u.indices[1:Lshifted]
            y_dir=rand([L+1, 1+L*L-L])            #chooses between [+a2 ,-a2]
            for iy in 1:L
                src=column+(iy-1)*L
                trg=Bsrcdir2trg[src, y_dir] # dir=L+1 is a shift ty=a*ey
                for slice in 1:nslices(mc)
                    c[:, trg, slice] .= tc[:, src, slice]
                end
            end
        end
    end
    return nothing
end




"""
    propose function for SpatialShuffle([mc, model])

        A global update that randomly swaps spatial indices of a configuration without 
        changing the temporal indices.
"""
@bm function propose_conf!(u::SpatialShuffle, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    shuffle!(u.indices)
    for slice in 1:nslices(mc), (i, j) in enumerate(u.indices)
        c[:,i, slice] .= tc[:, j, slice]
    end
    return nothing
end
"""
    propose function for TemporalShuffle([mc, model])

        A global update that randomly swaps time slice indices of a configuration 
        without changing the spatial indices.
"""
@bm function propose_conf!(u::TemporalShuffle, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    shuffle!(u.indices)
    for (k, l) in enumerate(u.indices), i in 1:length(lattice(mc))
        c[:, i, k] .= tc[:, i, l]
    end
    return nothing
end
"""
    SpatialStaggeredFlip([mc, model])

A global update that multiplies each spatial sheat with e^(i*(Q1+Q2)⋅r)=(-1)^(ix+iy).
"""
mutable struct SpatialStaggeredFlip <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    SpatialStaggeredFlip(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
SpatialStaggeredFlip(p, mc::DQMC, model::Model) = SpatialStaggeredFlip(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
SpatialStaggeredFlip(mc::DQMC, model::Model) = SpatialStaggeredFlip(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
SpatialStaggeredFlip() = SpatialStaggeredFlip(1.0,  0, 10000)

name(::SpatialStaggeredFlip) = "SpatialStaggeredFlip"
function _save(f::FileLike, name::String, update::SpatialStaggeredFlip)
    write(f, "$name/tag", :SpatialStaggeredFlip)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:SpatialStaggeredFlip}) = SpatialStaggeredFlip(f["p"], f["count"], f["full_update_interval"])

function init!(mc, u::SpatialStaggeredFlip)
    resize!(u.indices, length(lattice(mc)))
    u.indices .= 1:length(lattice(mc))
    nothing
end

@bm function propose_conf!(u::SpatialStaggeredFlip, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    L= model.l.Ls[1]
    u.count+=1;
    N=length(lattice(mc))
    k=Int8(length(field.x))   #number of discretization (typically k=4)
    if mod(u.count, u.full_update_interval)==0
        Nflipped=N
    else
        Nflipped=max(1, round(Int, u.p *N))
    end
    shuffle!(u.indices)
    for slice in 1:nslices(mc)
        for  i in u.indices[1:Nflipped]
            iy, ix=fldmod1(i, L) 
            if isodd(ix+iy) 
                c[:, i, slice] .= Int8(k+1) .- c[:, i, slice]
            end
        end
    end
    return nothing
end

"""
    SpatialStripedXorYFlip([mc, model])

A global update that multiplies each spatial sheet with e^(i*(Q1)⋅r)=(-1)^(ix).
"""
mutable struct SpatialStripedXorYFlip <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    SpatialStripedXorYFlip(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
SpatialStripedXorYFlip(p, mc::DQMC, model::Model) = SpatialStripedXorYFlip(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
SpatialStripedXorYFlip(mc::DQMC, model::Model) = SpatialStripedXorYFlip(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
SpatialStripedXorYFlip() = SpatialStripedXorYFlip(1.0,  0, 10000)

name(::SpatialStripedXorYFlip) = "SpatialStripedXorYFlip"
function _save(f::FileLike, name::String, update::SpatialStripedXorYFlip)
    write(f, "$name/tag", :SpatialStripedXorYFlip)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:SpatialStripedXorYFlip}) = SpatialStripedXorYFlip(f["p"], f["count"], f["full_update_interval"])

function init!(mc, u::SpatialStripedXorYFlip)
    resize!(u.indices, length(lattice(mc)))
    u.indices .= 1:length(lattice(mc))
    nothing
end

@bm function propose_conf!(u::SpatialStripedXorYFlip, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    L= model.l.Ls[1]
    u.count+=1;
    N=length(lattice(mc))
    k=Int8(length(field.x))   #number of discretization (typically k=4)

    if mod(u.count, u.full_update_interval)==0
        Nflipped=N
    else
        Nflipped=max(1, round(Int, u.p *N))
    end
    shuffle!(u.indices)
    if isodd(div(u.count, 3))   #trying three times x, then 3 times y, etc.
        for slice in 1:nslices(mc)
            for  i in u.indices[1:Nflipped]
                iy, ix=fldmod1(i, L) 
                if isodd(ix) 
                    c[:, i, slice] .= Int8(k+1) .- c[:, i, slice]
                end
            end
        end
    else
        for slice in 1:nslices(mc)
            for  i in u.indices[1:Nflipped]
                iy, ix=fldmod1(i, L) 
                if isodd(iy) 
                    c[:, i, slice] .= Int8(k+1) .- c[:, i, slice]
                end
            end
        end
    end
    return nothing
end

#### TODO: Only works for k=4!
"""
    computes the (probabilistic) outcome of (ϕ(i1) +ϕ(i2))/2
"""
function _add_field(i1::Int8, i2::Int8)
    i1==i2 && return i1
    ((i1==1 && i2==4) || (i1==4 && i2==1)) && return rand()<0.5 ? Int8(1) : Int8(4)
    isodd(i1-i2) && return rand()<0.5 ? i1 : i2 #adjacent values
    ((i1==1 && i2==3) || (i1==3 && i2==1)) && return rand()<0.5*sqrt(2+sqrt(3)) ? Int8(2) : Int8(1)
    ((i1==2 && i2==4) || (i1==4 && i2==2)) && return rand()<0.5*sqrt(2+sqrt(3)) ? Int8(3) : Int8(4)
    nothing
end

"""
    AddShiftedConfiguration([mc, model])
TODO: Only works for k=4!
A global update that tests for ϕ`(r)=[ϕ(r) + ϕ(r+tₓ)]/2, or ty.
Apparently, this update changes the weights [spin flip don`t do that],
which is problematic as it leads to fac ~ 10^4, and acc_rate=1. 
Possibly, violates detailed balance ?!
"""
mutable struct AddShiftedConfiguration <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    AddShiftedConfiguration(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
AddShiftedConfiguration(p, mc::DQMC, model::Model) = AddShiftedConfiguration(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
AddShiftedConfiguration(mc::DQMC, model::Model) = AddShiftedConfiguration(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
AddShiftedConfiguration() = AddShiftedConfiguration(1.0,  0, 10000)

name(::AddShiftedConfiguration) = "AddShiftedConfiguration"
function _save(f::FileLike, name::String, update::AddShiftedConfiguration)
    write(f, "$name/tag", :AddShiftedConfiguration)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:AddShiftedConfiguration}) = AddShiftedConfiguration(f["p"], f["count"], f["full_update_interval"])

function init!(mc, u::AddShiftedConfiguration)
    # resize!(u.indices, length(lattice(mc)))
    # u.indices .= 1:length(lattice(mc))
    L= mc.model.l.Ls[1]
    resize!(u.indices, L)
    u.indices .= 1:L
    nothing
end

@bm function propose_conf!(u::AddShiftedConfiguration, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    L= model.l.Ls[1]
    u.count+=1;
    dirx=isodd(div(u.count, 3))   #trying three times x, then 3 times y, etc.
    Bsrcdir2trg = model.l[:Bravais_srcdir2trg]::Matrix{Int}
    N=length(lattice(mc))
    # mod(u.count, u.full_update_interval)==0 ? Nchanged=N : Nchanged=max(1, round(Int, u.p *N))
    mod(u.count, u.full_update_interval)==0 ? Lchanged=L : Lchanged=max(1, round(Int, u.p *L))

    shuffle!(u.indices)
    for slice in 1:nslices(mc)
        # for  i in u.indices[1:Nchanged]
        #     dirx ? ip=Bsrcdir2trg[i,2] : ip=Bsrcdir2trg[i,L+1] 
        #     # dir=2 is a shift tx=a*ex
        #     # dir=5 is a shift ty=a*ey
        #     c[:, i, slice] .= _add_field.(tc[:, i, slice], tc[:, ip, slice])
        # end
        for  row in u.indices[1:Lchanged]
            for index in 1:L
                if dirx
                    i=index +(row-1)*L
                    ip=Bsrcdir2trg[i,2]
                else
                    i=row +(index-1)*L
                    ip=Bsrcdir2trg[i,L+1]
                end
                c[:, i, slice] .= _add_field.(tc[:, i, slice], tc[:, ip, slice])
            end
        end
    end
    return nothing
end

"""
    AddStaggeredConfiguration([mc, model])
TODO: Only works for k=4!

A global update that tests for ϕ`(r)=ϕ(r) *(1+e^(i(Q₁+Q₂)⋅r))/2.
"""
mutable struct AddStaggeredConfiguration <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    AddStaggeredConfiguration(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
AddStaggeredConfiguration(p, mc::DQMC, model::Model) = AddStaggeredConfiguration(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
AddStaggeredConfiguration(mc::DQMC, model::Model) = AddStaggeredConfiguration(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
AddStaggeredConfiguration() = AddStaggeredConfiguration(1.0,  0, 10000)

name(::AddStaggeredConfiguration) = "AddStaggeredConfiguration"
function _save(f::FileLike, name::String, update::AddStaggeredConfiguration)
    write(f, "$name/tag", :AddStaggeredConfiguration)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:AddStaggeredConfiguration}) = AddStaggeredConfiguration(f["p"], f["count"], f["full_update_interval"])

function init!(mc, u::AddStaggeredConfiguration)
    resize!(u.indices, length(lattice(mc)))
    u.indices .= 1:length(lattice(mc))
    nothing
end

@bm function propose_conf!(u::AddStaggeredConfiguration, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    L= model.l.Ls[1]
    u.count+=1;
    add=isodd(div(u.count, 3))   #trying three times x, then 3 times y, etc.
    N=length(lattice(mc))
    mod(u.count, u.full_update_interval)==0 ? Nchanged=N : Nchanged=max(1, round(Int, u.p *N))

    shuffle!(u.indices)
    for slice in 1:nslices(mc)
        for  i in u.indices[1:Nchanged]
            iy, ix=fldmod1(i, L)                          
            if add
                if isodd(ix+iy)
                    c[:, i, slice] .= _add_field.(tc[:, i, slice], (Int8(5) .- tc[:, i, slice]))
                end
            else
                if iseven(ix+iy)
                    c[:, i, slice] .= _add_field.(tc[:, i, slice], (Int8(5) .- tc[:, i, slice]))
                end
            end
        end
    end
    return nothing
end


"""
    LinWeightedStaggFlip([mc, model])

This global update is designed to promote the transition from the double-Q charge-density wave state
to one of the nematic states. In particular, it iterates over staggered lattice sites, 
flips them according to the two neighbors in x- or y-directions.
"""
mutable struct LinWeightedStaggFlip <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    LinWeightedStaggFlip(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
LinWeightedStaggFlip(p, mc::DQMC, model::Model) = LinWeightedStaggFlip(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
LinWeightedStaggFlip(mc::DQMC, model::Model) = LinWeightedStaggFlip(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
LinWeightedStaggFlip() = LinWeightedStaggFlip(1.0,  0, 10000)

name(::LinWeightedStaggFlip) = "LinWeightedStaggFlip"
function _save(f::FileLike, name::String, update::LinWeightedStaggFlip)
    write(f, "$name/tag", :LinWeightedStaggFlip)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:LinWeightedStaggFlip}) = LinWeightedStaggFlip(f["p"], f["count"], f["full_update_interval"])

function init!(mc, u::LinWeightedStaggFlip)
    resize!(u.indices, length(lattice(mc)))
    u.indices .= 1:length(lattice(mc))
    nothing
end

@bm function propose_conf!(u::LinWeightedStaggFlip, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    L= model.l.Ls[1]
    k=Int8(length(field.x))   #number of discretization (typically k=4)
    u.count+=1;
    odds=true
    dirx=isodd(div(u.count, 3))   #trying three times x, then 3 times y, etc.
    N=length(lattice(mc))
    Bsrcdir2trg = model.l[:Bravais_srcdir2trg]::Matrix{Int}

    mod(u.count, u.full_update_interval)==0 ? Nchanged=N : Nchanged=max(1, round(Int, u.p *N))

    shuffle!(u.indices)
    for slice in 1:nslices(mc)
        for  i in u.indices[1:Nchanged]
            iy, ix=fldmod1(i, L)    
            if dirx
                iP=Bsrcdir2trg[i, 2]     # i + a1
                iM=Bsrcdir2trg[i, L]     # i - a1   
            else
                iP=Bsrcdir2trg[i, 1 + L]     # i + a2
                iM=Bsrcdir2trg[i, 1+L*L-L]   # i - a2  
            end                      
            if (odds && isodd(ix+iy)) || (!(odds) && iseven(ix+iy))
                av_vec=0.5*(field.x[tc[:, iP, slice]]+ field.x[tc[:, iM, slice]])  
                if av_vec ⋅ field.x[tc[:, i, slice]] <0 
                    c[:, i, slice] .= Int8(k+1) .- tc[:, i, slice]                                              
                end
            end
        end
    end
    return nothing
end

"""
    LinWeightedFlip([mc, model])

This global update is designed to promote the transition from the double-Q charge-density wave state
to one of the nematic states. In particular, it iterates over randomly chosen lattice sites, 
and flips them according to the two neighbors in x- or y-directions.
"""
mutable struct LinWeightedFlip <: AbstractGlobalUpdate 
    p::Float64
    indices::Vector{Int}
    count::Int
    full_update_interval::Int
    LinWeightedFlip(p, count, full_update_interval) = new(p, Vector{Int}(undef, 0), count, full_update_interval)
end
LinWeightedFlip(p, mc::DQMC, model::Model) = LinWeightedFlip(p, 0, mc.parameters.sweeps+mc.parameters.thermalization)
LinWeightedFlip(mc::DQMC, model::Model) = LinWeightedFlip(1.0,  0, mc.parameters.sweeps+mc.parameters.thermalization)
LinWeightedFlip() = LinWeightedFlip(1.0,  0, 10000)

name(::LinWeightedFlip) = "LinWeightedFlip"
function _save(f::FileLike, name::String, update::LinWeightedFlip)
    write(f, "$name/tag", :LinWeightedFlip)
    write(f, "$name/p", update.p)
    write(f, "$name/count", update.count)
    write(f, "$name/full_update_interval", update.full_update_interval)
end
_load(f::FileLike, ::Val{:LinWeightedFlip}) = LinWeightedFlip(f["p"], f["count"], f["full_update_interval"])


function init!(mc, u::LinWeightedFlip)
    resize!(u.indices, length(lattice(mc)))
    u.indices .= 1:length(lattice(mc))
    nothing
end

@bm function propose_conf!(u::LinWeightedFlip, mc, model, field::AbstractDiscreteMBF)
    c = conf(field); tc = temp_conf(field)
    copyto!(tc, c)
    L= model.l.Ls[1]
    k=Int8(length(field.x))   #number of discretization (typically k=4)
    u.count+=1;
    dirx=isodd(div(u.count, 3))   #trying three times x, then 3 times y, etc.
    N=length(lattice(mc))
    Bsrcdir2trg = model.l[:Bravais_srcdir2trg]::Matrix{Int}

    mod(u.count, u.full_update_interval)==0 ? Nchanged=N : Nchanged=max(1, round(Int, u.p *N))

    shuffle!(u.indices)
    for slice in 1:nslices(mc)
        for  i in u.indices[1:Nchanged]
            if dirx
                iP=Bsrcdir2trg[i, 2]     # i + a1
                iM=Bsrcdir2trg[i, L]     # i - a1   
            else
                iP=Bsrcdir2trg[i, 1 + L]     # i + a2
                iM=Bsrcdir2trg[i, 1+L*L-L]   # i - a2  
            end                      
            av_vec=0.5*(field.x[tc[:, iP, slice]]+ field.x[tc[:, iM, slice]])  
            if av_vec ⋅ field.x[tc[:, i, slice]] <0 
                c[:, i, slice] .= Int8(k+1) .- tc[:, i, slice]                                              
            end 
        end
    end
    return nothing
end