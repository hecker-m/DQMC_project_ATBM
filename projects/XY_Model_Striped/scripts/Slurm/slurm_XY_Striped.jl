#!/usr/bin/bash -l
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5g
#SBATCH --mail-type=all
#SBATCH --mail-user=mhecker@umn.edu
#SBATCH --array=1-160
#SBATCH --job-name=XY_Striped
#SBATCH -o %x-%A_%a.out
#=
    pwd
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    echo
    srun julia --threads=$SLURM_CPUS_PER_TASK slurm_XY_Striped.jl
    exit
=#

include("../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(1)

@show jobname = ENV["SLURM_JOB_NAME"]
@show jobid = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])


# Note that we implement U=U/Nϕ (see below)

Us=[0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.3, 1.6, 2.0, 2.4]
βs=[1, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 10, 20]
paras=[(L=8, β=β0, U=U0, Pe=true) for U0 in Us, β0 in βs][:]

therm = 1000
sweeps = 2^(13)
L=paras[jobid].L;
beta=paras[jobid].β;
U=paras[jobid].U;
peierls=paras[jobid].Pe;
N=L^2;
T=1/beta

t=1.0;
t1a=0.6*t;
t1b=0.2*t;
tpa=0.8*t;
tpb=1.0*t;
μ0=1.0*t;
ϵ0=0.0;
ϵB1Plus=ϵ0*t1a;
ϵB1Minus=ϵ0*t1b;



@show nworkers=parse(Int, ENV["SLURM_CPUS_PER_TASK"])

path1=pwd() * "/run_saves/eps_$(Int(100ϵ0))pc";
path=path1 * "/D_L$(L)_b_" * to_string(beta) * "_U_" * to_string(U) * "_eps_$(Int(100ϵ0))pc_B_$(Int(peierls))";

sleep(2jobid)
if !isdir(path1)
    mkdir(path1)
end
if !isdir(path)
    mkdir(path)
end


@inbounds Threads.@threads for worker in 1:nworkers

    #######
    #Choosing the Configuration field
    #######
    Nϕ=2            #The code is implemented to work for Nϕ=1,2,3, i.e. Ising, XY-, or Heisenberg field
    discrete=true  #true-> discrete field, false->continuous field
    if discrete
        prefix="Discrete"
        if Nϕ==1
            suffix="_X_symm";
        elseif Nϕ==2
            suffix="_symm";
        else
            suffix=""; #by default we use the symmetry-optimized Ising- and XY- fields
        end
        schedule=SimpleScheduler(LocalSweep(),PartialGlobalFlip(0.2, 0, 10), LocalSweep(8),
            LocalSweep(),SpatialStaggeredFlip(0.2, 2, 10), LocalSweep(8),
            LocalSweep(),pGlobalXorYshift(0.2, 4, 10), LocalSweep(8),
            LocalSweep(),AddStaggeredConfiguration(0.2, 6, 10), LocalSweep(8),
            LocalSweep(),LinWeightedFlip(0.2, 8, 10), LocalSweep(8))
    else
        prefix="Cont";suffix="";
        schedule=SimpleScheduler(GlobalConstMove(1), LocalSweep(5))
    end
    field=Symbol(prefix * "_MBF$(Nϕ)" *suffix)

    #######
    #Initializing the Model, and the DQMC Simulation
    #######
    model = TwoBandModel(dims=2, L=L, U = U/Nϕ, tx = [t1a+ϵB1Plus; t1b+ϵB1Minus], ty = [t1a-ϵB1Plus; t1b-ϵB1Minus],
        tp = [tpa; tpb], μs = [μ0; -μ0], peierls=peierls)


    mc = DQMC(model, field=getfield(Main, field), Nϕ=Nϕ, beta=beta, scheduler = schedule,
            safe_mult = 5, delta_tau=0.05, measure_rate = 10, box_local = 15.0, box_global = 1.0,
            recording_rate=10)

    #######
    #Adding thermalization measurements
    #######
    metype=ComplexF64
    #mc.thermalization_measurements[:occ_th] = occupation(mc, model,obs=FullBinner(Vector{metype}))
    #mc.thermalization_measurements[:E_th] = total_energy(mc, model, obs=FullBinner(metype))
    #mc.thermalization_measurements[:Mx_z_th] = magnetization(mc, model, :z, kernel = Mx_z_kernel, obs=FullBinner(Vector{metype}))



    ##############
    # Adding measurements
    ##############
    cap=16*ceil(Int, sweeps /mc.parameters.measure_rate);
    mc.measurements[:occ] = occupation(mc, model, capacity=cap)
    mc.measurements[:E] = total_energy(mc, model, obs=FullBinner(Float64))

    ##### magnetization
    #
    mc.measurements[:Mx_X_OP] = Mx_X_OP(mc, model, :x)
    mc.measurements[:Mx_Y_OP] = Mx_X_OP(mc, model, :y)
    mc.measurements[:Mx_Z_OP] = Mx_X_OP(mc, model, :z)

    mc.measurements[:SDC_Mx_z] = spin_density_correlation(mc, model, :z, kernel=full_sdc_Mx_z_kernel, capacity=cap)
    mc.measurements[:SDS_Mx_z] = spin_density_susceptibility(mc, model, :z, kernel=full_sdc_Mx_z_kernel, capacity=cap)
    mc.measurements[:SDC_Mx_x] = spin_density_correlation(mc, model, :x, kernel=full_sdc_Mx_x_kernel, capacity=cap)
    mc.measurements[:SDS_Mx_x] = spin_density_susceptibility(mc, model, :x, kernel=full_sdc_Mx_x_kernel, capacity=cap)
    mc.measurements[:SDC_Mx_y] = spin_density_correlation(mc, model, :y, kernel=full_sdc_Mx_y_kernel, capacity=cap)
    mc.measurements[:SDS_Mx_y] = spin_density_susceptibility(mc, model, :y, kernel=full_sdc_Mx_y_kernel, capacity=cap)

    ##### superconductivity
    #
    #mc.measurements[:PDC_s] = pairing_correlation(mc, model, kernel = pc_swave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_s] = pairing_susceptibility(mc, model, kernel = pc_swave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    #mc.measurements[:PDC_spm] = pairing_correlation(mc, model, kernel = pc_spm_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_spm] = pairing_susceptibility(mc, model, kernel = pc_spm_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_XX] = pairing_susceptibility(mc, model, kernel = pc_XX_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_YYzz] = pairing_susceptibility(mc, model, kernel = pc_YYzz_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)

    mc.measurements[:Δ_Zy_bil_OP] = Δ_Zy_bil_OP(mc, model)
    mc.measurements[:Δ_0y_bil_OP] = Δ_0y_bil_OP(mc, model)
    mc.measurements[:Δ_Xy_bil_OP] = Δ_Xy_bil_OP(mc, model)
    mc.measurements[:Δ_Ysum_bil_OP] = Δ_Ysum_bil_OP(mc, model)


    ##### charge sector
    #
    #mc.measurements[:CDC] = charge_density_correlation(mc, model, capacity=cap)
    mc.measurements[:CDS] = charge_density_susceptibility(mc, model, capacity=cap)
    #mc.measurements[:Zk_proxy] = spectral_weight_proxy(mc, model, capacity=cap)
    #mc.measurements[:green]= greens_measurement(mc, model, capacity=cap)
    #mc.measurements[:CD_X_OP] = CD_X_OP(mc, model)  #For reasons unclear, the ρ^{x0} -OP is purely imaginary
    #CD_X OP is zero in case of Ising-X and XY
    mc.measurements[:CDSxx]=charge_density_susceptibility(mc, model, capacity=cap, flavor_iterator = MonteCarlo.FlavorIterator(mc, 0),
            kernel = full_cdc_XX_kernel)


    ##### B1-nematic sector
    #
    mc.measurements[:B1_OP] = nematic_OP(mc, model, obs=FullBinner(Float64))
    mc.measurements[:B1_proxy_OP] = proxy_B1_OP(mc, model)
    mc.measurements[:B1_CDS] = B1_charge_density_susceptibility(mc, model, capacity=cap)
    mc.measurements[:NemS] = nematic_susceptibility(mc, model, capacity=cap)
    mc.measurements[:NemC] = nematic_correlation(mc, model, capacity=cap)

    ##### A1p double-Q sector
    #
    mc.measurements[:A1p_OP] = A1p_OP(mc, model, obs=FullBinner(Float64))
    mc.measurements[:A1p_proxy_OP] = proxy_A1p_OP(mc, model)
    mc.measurements[:A1p_dQ_S] = A1_Q1Q2_susceptibility(mc, model, capacity=cap)
    mc.measurements[:A1p_dQ_C] = A1_Q1Q2_correlation(mc, model, capacity=cap)

    ##### B1p double-Q sector
    #
    #mc.measurements[:B1p_OP_x] = B1p_OP(mc, model, :x, obs=FullBinner(Float64))
    #mc.measurements[:B1p_OP_y] = B1p_OP(mc, model, :y, obs=FullBinner(Float64))
    mc.measurements[:B1p_OP_z] = B1p_OP(mc, model, :z, obs=FullBinner(Float64))
    #B1p OP is zero in x, and y direction in case of Ising-X
    mc.measurements[:B1p_dQ_C] = B1p_Q1Q2_correlation(mc, model, capacity=cap)
    mc.measurements[:B1p_dQ_S] = B1p_Q1Q2_susceptibility(mc, model, capacity=cap)


    ################
    # run the simulation
    ###############

    st="D_XY_b_" * to_string(beta) *"_U_"* to_string(U) *"_L$(L)_eps_$(Int(100ϵ0))pc_B_$(Int(peierls))_sw$(sweeps)_th$(therm)_worker_$(worker).jld2";
    resumable_file=path* "/resumable_" *st;
    final_file=path * "/" * st;

    exit_code= run!(mc, sweeps=sweeps, thermalization=therm, verbose=false, safe_every= Hour(2),  ten_reg=3,
            resumable_filename= resumable_file, overwrite=true)

    MonteCarlo.save(final_file, mc, overwrite = true, rename = false)

    if Int(exit_code)==0 && @isdefined(resumable_file) && isfile(resumable_file)
        rm(resumable_file)
    end

end


