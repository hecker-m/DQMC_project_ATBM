#!/usr/bin/bash -l
#SBATCH --time=28:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5g
#SBATCH --mail-type=all
#SBATCH --mail-user=mhecker@umn.edu
#SBATCH --array=1-6
#SBATCH --job-name=IsingX_Neel
#SBATCH -o %x-%A_%a.out
#=
    pwd
    echo $SLURM_NPROCS
    echo $SLURM_CPUS_PER_TASK
    echo
    srun julia --threads=$SLURM_CPUS_PER_TASK slurm_IsingX_Neel.jl
    exit
=#

include("../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(1)

@show jobname = ENV["SLURM_JOB_NAME"]
@show jobid = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])


# Note that we implement U=U/Nϕ (see below)

Us=[0.6, 0.7, 0.8, 1.0, 1.2, 1.3]
βs=[10, ]
paras=[(L=8, β=β0, U=U0, Pe=true) for U0 in Us, β0 in βs][:]

therm = 1000
sweeps = 2^(13)
L=paras[jobid].L;
beta=paras[jobid].β;
U=paras[jobid].U;
peierls=paras[jobid].Pe;
N=L^2;
T=1/beta

t=1.0;δ=0.4*t;μ0=-2*t;

@show nworkers=parse(Int, ENV["SLURM_CPUS_PER_TASK"])

path=pwd() * "/run_saves/D_L$(L)_b_" * to_string(beta) * "_U_" * to_string(U) * "_B_$(Int(peierls))";
sleep(2jobid)
if !isdir(path)
    mkdir(path)
end


@inbounds Threads.@threads for worker in 1:nworkers

    #######
    #Choosing the Configuration field
    #######
    Nϕ=1            #The code is implemented to work for Nϕ=1,2,3, i.e. Ising, XY-, or Heisenberg field
    discrete=true  #true-> discrete field, false->continuous field
    if discrete
        prefix="Discrete"
        if Nϕ==1
            suffix="_X_symm";
        else
            suffix=""; #by default we use the symmetry-optimized Ising field
        end
        schedule=SimpleScheduler(LocalSweep())
    else
        prefix="Cont";suffix="";
        schedule=SimpleScheduler(GlobalConstMove(1), LocalSweep(5))
    end
    field=Symbol(prefix * "_MBF$(Nϕ)" *suffix)

    #######
    #Initializing the Model, and the DQMC Simulation
    #######
    model = TwoBandModel(dims=2, L=L, U = U/Nϕ, tx = [-(t+δ);-(t-δ)], ty = [-(t-δ);-(t+δ)],
        tp = [0.0; 0.0], μs = [μ0;-μ0], peierls=peierls)

    mc = DQMC(model, field=getfield(Main, field), Nϕ=Nϕ, beta=beta, scheduler = schedule,
            safe_mult = 5, delta_tau=0.05, measure_rate = 10, box_local = 15.0, box_global = 1.0,
            recording_rate=20)

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

    #mc.measurements[:SDC_Mx_z] = spin_density_correlation(mc, model, :z, kernel=full_sdc_Mx_z_kernel, capacity=cap)
    mc.measurements[:SDS_Mx_z] = spin_density_susceptibility(mc, model, :z, kernel=full_sdc_Mx_z_kernel, capacity=cap)
    #mc.measurements[:SDC_Mx_x] = spin_density_correlation(mc, model, :x, kernel=full_sdc_Mx_x_kernel, capacity=cap)
    mc.measurements[:SDS_Mx_x] = spin_density_susceptibility(mc, model, :x, kernel=full_sdc_Mx_x_kernel, capacity=cap)
    #mc.measurements[:SDC_Mx_y] = spin_density_correlation(mc, model, :y, kernel=full_sdc_Mx_y_kernel, capacity=cap)
    mc.measurements[:SDS_Mx_y] = spin_density_susceptibility(mc, model, :y, kernel=full_sdc_Mx_y_kernel, capacity=cap)

    ##### superconductivity
    #
    #mc.measurements[:PDC_s] = pairing_correlation(mc, model, kernel = pc_swave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_s] = pairing_susceptibility(mc, model, kernel = pc_swave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    #mc.measurements[:PDC_spm] = pairing_correlation(mc, model, kernel = pc_spm_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_spm] = pairing_susceptibility(mc, model, kernel = pc_spm_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)

    ##### charge sector
    #
    #mc.measurements[:CDC] = charge_density_correlation(mc, model, capacity=cap)
    mc.measurements[:CDS] = charge_density_susceptibility(mc, model, capacity=cap)
    #mc.measurements[:Zk_proxy] = spectral_weight_proxy(mc, model, capacity=cap)


    ################
    # run the simulation
    ###############

    st="DSx_FP_b_" * to_string(beta) *"_U_"* to_string(U) *"_L$(L)_B_$(Int(peierls))_sw$(sweeps)_th$(therm)_worker_$(worker).jld2";
    resumable_file=path* "/resumable_" *st;
    final_file=path * "/" * st;

    exit_code= run!(mc, sweeps=sweeps, thermalization=therm, verbose=false, safe_every= Hour(2),  ten_reg=3,
            resumable_filename= resumable_file, overwrite=true)

    MonteCarlo.save(final_file, mc, overwrite = true, rename = false)

    if Int(exit_code)==0 && @isdefined(resumable_file) && isfile(resumable_file)
        rm(resumable_file)
    end

end


