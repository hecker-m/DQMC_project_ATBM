module MonteCarlo

using Octavian
using LinearAlgebra: AbstractMatrix, Hermitian
using Reexport
# Loading the RNG will fail if Random is nto exported
@reexport using Random, BinningAnalysis
using Parameters, Requires
using TimerOutputs, StructArrays
using Printf, SparseArrays, LinearAlgebra, Dates, Statistics, Random, Distributed
using FFTW
# unoptized version, mainly for tests
# this still uses custom linear algebra methods, just without loopvectorization
if get(ENV, "MONTECARLO_USE_LOOPVECTORIZATION", "true") == "true"
    import LoopVectorization
    using LoopVectorization: @turbo
else
    printstyled(
        "Using MonteCarlo.jl without LoopVectorization. This should only be done for tests.",
        color = :red
    )
    macro turbo(code)
        esc(quote @inbounds @fastmath $code end)
    end
end

import JLD2, CodecLz4

include("helpers.jl")
export enable_benchmarks, disable_benchmarks, print_timer, reset_timer!

mutable struct FileData
    data::Dict{String, Any}
    subspace::String
    path::String
end
FileData(data, path) = FileData(data, "", path)
@bm function Base.getindex(f::FileData, key::String)
    subspace = f.subspace * key * '/' 
    if any(k -> startswith(k, subspace), keys(f.data))
        return FileData(f.data, subspace, f.path)
    elseif haskey(f.data, f.subspace * key)
        return getindex(f.data, f.subspace * key)
    else
        return getindex(f.data, key)
    end
end
@bm function Base.haskey(f::FileData, key::String)
    subspace = f.subspace * key
    haskey(f.data, key) || any(k -> startswith(k, subspace), keys(f.data))
end
@bm function Base.get(f::FileData, k::String, v)
    return haskey(f, k) ? getindex(f, k) : v
end
@bm function Base.keys(f::FileData)
    _keys = String[]
    i0 = length(f.subspace)
    for key in keys(f.data)
        if startswith(key, f.subspace)
            i1 = something(findnext(isequal('/'), key, i0+1), length(key)+1)
            _key = string(key[i0+1:i1-1])
            _key in _keys || push!(_keys, _key)
        end
    end
    return _keys
end

const FileLike = Union{JLD2.JLDFile, JLD2.Group, FileData}
filepath(f::FileLike) = f.path
filepath(g::JLD2.Group) = g.f.path


include("lattices/lattice.jl")
include("flavors/abstract.jl")
include("models/abstract.jl")

include("configurations.jl")
export Discarder, ConfigRecorder, BufferedConfigRecorder, RelativePath, AbsolutePath
include("Measurements.jl")
export measurements, observables

include("lattices/constructors.jl")
include("lattices/lattice_cache.jl")
include("lattices/lattice_iterators.jl")
include("lattices/directions.jl")

export AbstractLattice, Lattice, Bravais
export Chain, SquareLattice, CubicLattice, TriangularLattice, Honeycomb
# export AbstractLattice, Chain, SquareLattice, CubicLattice, TriangularLattice, ALPSLattice
export EachSite, EachSiteAndFlavor, OnSite, EachSitePair, EachSitePairByDistance, 
        EachLocalQuadByDistance, EachLocalQuadBySyncedDistance, 
        EachBondPairByBravaisDistance, Sum, ApplySymmetries
        
export unitcell, positions, bonds, bonds_open, directions

include("flavors/MC/MC.jl")
include("flavors/DQMC/main.jl")
export Greens, GreensAt, CombinedGreensIterator, TimeIntegral
export GreensMatrix, swapop, Restructure
export boson_energy_measurement, greens_measurement, occupation, occupation_summed, magnetization
export charge_density, charge_density_correlation, charge_density_susceptibility
export spin_density, spin_density_correlation, spin_density_susceptibility
export pairing, pairing_correlation, pairing_susceptibility
export current_current_susceptibility, superfluid_density
export kinetic_energy, interaction_energy, total_energy
export add_default_measurements!

export EmptyScheduler, SimpleScheduler, AdaptiveScheduler
export Adaptive, NoUpdate, LocalSweep
export GlobalFlip, GlobalShuffle, SpatialShuffle, TemporalShuffle, 
        Denoise, DenoiseFlip, StaggeredDenoise
export ReplicaExchange, ReplicaPull, connect, disconnect
export ChemicalPotentialTuning
# export mask, uniform_fourier, structure_factor, SymmetryWrapped, swave, eswave

include("models/Ising/IsingModel.jl")
include("models/HubbardModel.jl")
include("models/DummyModel.jl")
export IsingEnergyMeasurement, IsingMagnetizationMeasurement

include("FileIO.jl")
include("BinningAnalysis.jl")
export save, load, resume!


export reset!
export run!, resume!, replay!
export Model, IsingModel
# export RepulsiveGHQHubbardModel, AttractiveGHQHubbardModel
export HubbardModel, HubbardModelAttractive, HubbardModelRepulsive, AttractiveHubbardModel, RepulsiveHubbardModel
export MonteCarloFlavor, MC, DQMC, DQMCMeasurement
export greens, greens!, lattice, model, parameters
export DensityHirschField, MagneticHirschField, DensityGHQField, MagneticGHQField

# For extending
export AbstractMeasurement, Model


function __init__()
    # @require LatPhysBase="eec5c15a-e8bd-11e8-0d23-6799ca40c963" include("lattices/LatPhys.jl")
    @require LatticePhysics = "53011200-ee7a-11e8-39f1-5f3e57afe4fd" include("lattices/LatPhys.jl")
    @require MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195" begin
        include("mpi.jl")
        export mpi_queue
        
        include("flavors/DQMC/updates/mpi_updates.jl")
        export MPIReplicaExchange, MPIReplicaPull
    end

    return nothing
end
include("my_files/my_linalg.jl")
include("my_files/helpful_fcns.jl")
include("my_files/TwoBandModel_import.jl")
include("my_files/TBM_field_import.jl")
include("my_files/TBM_field_discrete.jl")
include("my_files/TBM_symm_field.jl")
include("my_files/my_lattice_iterators.jl")

include("my_files/measurements/spectral_weight_proxy.jl")
include("my_files/my_kernels.jl")
include("my_files/measurements/nematic_susceptibility.jl")
include("my_files/measurements/doubleQ_A1p_susceptibility.jl")
include("my_files/measurements/doubleQ_B1p_susceptibility.jl")

include("my_files/measurements/B1_charge_density_susceptibility.jl")
include("my_files/measurements/order_parameters.jl")
include("my_files/measurements/diamagnetic_Kx.jl")
include("my_files/measurements/PhaseStiffness_measurement.jl")
include("my_files/measurements/heat_capacity.jl")
include("my_files/q_discretization.jl")
include("my_files/my_updates.jl")

#utility functions
export _zero!, to_string

#updates that I defined
export GlobalConstMove, SpatialStaggeredFlip, SpatialStripedXorYFlip, PartialGlobalFlip, pGlobalXorYshift, AddShiftedConfiguration, AddStaggeredConfiguration, LinWeightedStaggFlip, LinWeightedFlip
# Fields and model that I defined
export TwoBandModel, hopping_matrix, AbstractMagnBosonField, AbstractContMBF, Cont_MBF1, Cont_MBF2, Cont_MBF3, AbstractDiscreteMBF, Discrete_MBF1, Discrete_MBF1_symm, Discrete_MBF1_X, Discrete_MBF1_X_symm,  Discrete_MBF2, Discrete_MBF2_symm, Discrete_MBF3

#measurements that I defined
export nematic_measurement, nematic_correlation, nematic_susceptibility, full_nem_kernel
export A1_Q1Q2_measurement, A1_Q1Q2_correlation, A1_Q1Q2_susceptibility, full_A1_Q1Q2_kernel
export B1_charge_density_measurement, B1_charge_density_correlation, B1_charge_density_susceptibility, full_B1_charge_density_kernel
export nematic_OP_measurement, nematic_OP, full_nematic_OP_kernel, A1p_OP_measurement, A1p_OP, full_B1p_OP_kernel, B1p_OP, B1p_OP_measurement
export Δ_Zy_bil_OP, Δ_0y_bil_OP, Δ_Ysum_bil_OP, Δ_Xy_bil_OP
export proxy_A1p_OP, proxy_A1p_OP_kernel, proxy_B1_OP, proxy_B1_OP_kernel, CD_X_OP, charge_X_OP_kernel
export B1p_Q1Q2_measurement, full_B1p_Q1Q2_kernel, B1p_Q1Q2_correlation, B1p_Q1Q2_susceptibility

export kinetic_energy_kernel, interaction_energy_kernel, total_energy_kernel, Mx_z_kernel, Mx_x_kernel, Mx_X_OP, full_sdc_Mx_z_kernel, full_sdc_Mx_x_kernel, full_sdc_Mx_y_kernel, full_cdc_kernel, reduced_cdc_kernel, reduced_cdc_kernel_old, pc_combined_kernel, nearest_neighbor_count, pc_swave_kernel, pc_swave_kernel_conj, pc_swave_kernel_symm, pc_spm_wave_kernel, pc_spm_wave_kernel_conj, pc_spm_wave_kernel_symm, pc_XX_wave_kernel, pc_YYzz_wave_kernel, my_cc_kernel, _reliable_level, spectral_weight_proxy, full_sdc_Mx_z_kernel2, full_cdc_XX_kernel
export kx_dia_measurement, phase_stiffness
export heat_cap_h2, heat_cap_h3, heat_cap_h4, heat_cap_h2_kernel, heat_cap_h3_kernel, heat_cap_h4_kernel



export kinetic_energy_kernel_bonds
export FieldCache, unique_flavors, DQMCStack, init!, conf, nslices, interaction_matrix_exp!, interaction_matrix_exp_op!, calculate_detratio!, update_greens!, propose_local, randuniform, initialize_run

export _GM, FileLike, matrix_type, DQMCParameters, StandardFieldCache, FVec64, FMat64, AbstractField, field 
export interaction_matrix_type, hopping_matrix_type, greens_matrix_type, greens_eltype, interaction_eltype


#Iterators which I defined
export EachSitePair_summed, EachDoubleSitePairByDistance, EachDoubleSitePairByDistance_Q1Q2, EachSitePair_B1, EachSitePair_B1_OP, EachSitePair_A1p_OP, EachSitePair_B1p_OP, EachDoubleSitePairByDistance_B1p_Q1Q2, EachWeightedBond, PS_EachBondPairByBravaisDistance, EachDistancedBondPairSummed, EachBondEachSiteSummed, EachSiteTwiceSummed

#include("precompile.jl")
#include("Google Drive/DQMC/AFM_2_band_model/2BM_DQMC_code/TBM_field.jl")
#include("../../../TBM_field.jl")





end # module
