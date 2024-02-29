
# include("../../../../src/MonteCarlo.jl_modified/src/MonteCarlo2.jl")
using .MonteCarlo
using Distributions, DataFrames, JLD2, Dates, LinearAlgebra


# Note that we implement U=U/Nϕ (see below)

Us=[0.8, ]
βs=[1, ]
paras=[(L=4, β=β0, U=U0, Pe=true) for U0 in Us, β0 in βs][:]

jobid=1;
therm = 2000
sweeps = 2^(13)
L=paras[jobid].L;
beta=paras[jobid].β;
U=paras[jobid].U;
peierls=paras[jobid].Pe;
N=L^2;
T=1/beta

δτ=0.05;
t=0.0;
t0=1.0;
t1a=0.6*t;
t1b=0.2*t;
tpa=0.8*t;
tpb=1.0*t;
μ0=2.0 *t0;
ϵ0=0.0;
ϵB1Plus=ϵ0*t1a;
ϵB1Minus=ϵ0*t1b;
 
    #######
    #Choosing the Configuration field
    #######
    Nϕ=1            #The code is implemented to work for Nϕ=1,2,3, i.e. Ising, XY-, or Heisenberg field
    discrete=true  #true-> discrete field, false->continuous field
    k=8;
    if discrete
        k_string= k==4 ? "" : "_$(k)"
        prefix="Discrete" * k_string
        if Nϕ==1
            suffix="_X_symm";
        else
            suffix=""; #by default we use the symmetry-optimized Ising field
        end
        # schedule=SimpleScheduler(LocalSweep(),PartialGlobalFlip(0.2, 0, 10), LocalSweep(8),
        #     LocalSweep(),SpatialStaggeredFlip(0.2, 2, 10), LocalSweep(8),
        #     LocalSweep(),pGlobalXorYshift(0.2, 4, 10), LocalSweep(8),
        #     LocalSweep(),AddStaggeredConfiguration(0.2, 6, 10), LocalSweep(8),
        #     LocalSweep(),LinWeightedFlip(0.2, 8, 10), LocalSweep(8))
        schedule=SimpleScheduler(LocalSweep(),PartialGlobalFlip(0.2, 0, 10), LocalSweep(8),
            LocalSweep(),SpatialStaggeredFlip(0.2, 2, 10), LocalSweep(8),
            LocalSweep(),pGlobalXorYshift(0.2, 4, 10), LocalSweep(8))
        # schedule=SimpleScheduler(LocalSweep())
    else
        prefix="Cont";
        suffix="_X_symm";
        #schedule=SimpleScheduler(LocalSweep(), GlobalConstMove(1), LocalSweep(8))
        #schedule=SimpleScheduler(LocalSweep(), GlobalConstMove(1), LocalSweep(18))
        schedule=SimpleScheduler(LocalSweep())
    end
    #suffix="_symm";
    field=Symbol(prefix * "_MBF$(Nϕ)" *suffix)

    #######
    #Initializing the Model, and the DQMC Simulation
    #######
    model = TwoBandModel(dims=2, L=L, U = U/Nϕ, tx = [t1a+ϵB1Plus; t1b+ϵB1Minus], ty = [t1a-ϵB1Plus; t1b-ϵB1Minus],
        tp = [tpa; tpb], μs = [μ0; -μ0], peierls=peierls)


    mc = DQMC(model, field=getfield(Main, field), Nϕ=Nϕ, beta=beta, scheduler = schedule,
            safe_mult = 5, delta_tau=δτ, measure_rate = 10, box_local = 15.0, box_global = 1.0,
            recording_rate=10)

    #######
    #Adding thermalization measurements
    #######

    En_eltype=Float64;
    mc.thermalization_measurements[:E] = total_energy(mc, model, obs=FullBinner(En_eltype))
    mc.thermalization_measurements[:E_TI] = total_energy(mc, model, greens_iterator = TimeIntegral(mc),
         lattice_iterator = energy_iterator(), obs=FullBinner(En_eltype))

    mc.thermalization_measurements[:E_kin] = kinetic_energy(mc, model, obs=FullBinner(En_eltype))
    mc.thermalization_measurements[:E_pot] = interaction_energy(mc, model, obs=FullBinner(En_eltype))
        heat_eltype=Float64;
    
    mc.thermalization_measurements[:heat_cap_h2] = heat_cap_h2(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.thermalization_measurements[:heat_cap_h3] = heat_cap_h3(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.thermalization_measurements[:heat_cap_h4] = heat_cap_h4(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)

    mc.thermalization_measurements[:heat_cap_h2_TI] = heat_cap_h2_TI(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.thermalization_measurements[:heat_cap_h3_TI] = heat_cap_h3_TI(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.thermalization_measurements[:heat_cap_h4_TI] = heat_cap_h4_TI(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)

    mc.thermalization_measurements[:heat_cap_h4_onsite] = heat_cap_h4(mc, model, lattice_iterator=EachSiteSummed(), obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.thermalization_measurements[:heat_cap_h4_onsite_TI] = heat_cap_h4_TI(mc, model, lattice_iterator=EachSiteSummed(), obs=FullBinner(heat_eltype), eltype=heat_eltype)

    ##### magnetization
    mc.thermalization_measurements[:Mx_X_OP] = Mx_X_OP(mc, model, :x)
    mc.thermalization_measurements[:Mx_Y_OP] = Mx_X_OP(mc, model, :y)
    mc.thermalization_measurements[:Mx_Z_OP] = Mx_X_OP(mc, model, :z)

    mc.thermalization_measurements[:SDC_Mx_z] = spin_density_correlation(mc, model, :z, kernel=full_sdc_Mx_z_kernel, obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:SDS_Mx_z] = spin_density_susceptibility(mc, model, :z, kernel=full_sdc_Mx_z_kernel, obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:SDC_Mx_x] = spin_density_correlation(mc, model, :x, kernel=full_sdc_Mx_x_kernel, obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:SDS_Mx_x] = spin_density_susceptibility(mc, model, :x, kernel=full_sdc_Mx_x_kernel, obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:SDC_Mx_y] = spin_density_correlation(mc, model, :y, kernel=full_sdc_Mx_y_kernel, obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:SDS_Mx_y] = spin_density_susceptibility(mc, model, :y, kernel=full_sdc_Mx_y_kernel, obs=FullBinner(Array{ComplexF64}))

    mc.thermalization_measurements[:SDC_A1_Mx_x] = spin_density_correlation(mc, model, :x, kernel=full_sdc_Mx_x_kernel, 
        lattice_iterator = EachSitePair_A1(), obs=FullBinner(Float64))
    mc.thermalization_measurements[:SDS_A1_Mx_x] = spin_density_susceptibility(mc, model, :x, kernel=full_sdc_Mx_x_kernel, 
        lattice_iterator = EachSitePair_A1(), obs=FullBinner(Float64))

    ##### superconductivity
    mc.thermalization_measurements[:PDC_s] = pairing_correlation(mc, model, kernel = pc_swave_kernel, 
        lattice_iterator = EachSitePairByDistance(), obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:PDS_s] = pairing_susceptibility(mc, model, kernel = pc_swave_kernel, 
        lattice_iterator = EachSitePairByDistance(), obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:PDC_spm] = pairing_correlation(mc, model, kernel = pc_spm_wave_kernel, 
        lattice_iterator = EachSitePairByDistance(), obs=FullBinner(Array{ComplexF64}))
    mc.thermalization_measurements[:PDS_spm] = pairing_susceptibility(mc, model, kernel = pc_spm_wave_kernel, 
        lattice_iterator = EachSitePairByDistance(), obs=FullBinner(Array{ComplexF64}))

    # ##### charge sector
    mc.thermalization_measurements[:occ_summed] =occupation_summed(mc, model, eltype=Float64)
    mc.thermalization_measurements[:CDC_summed] = charge_density_correlation(mc, model, lattice_iterator = EachSitePair_summed(), 
        eltype =Float64, obs=FullBinner(Float64))

    # mc.thermalization_measurements[:CDS] = charge_density_susceptibility(mc, model, obs=FullBinner(Float64))
    mc.thermalization_measurements[:rho_s]  = phase_stiffness(mc, model, obs=FullBinner(Vector{ComplexF64}), eltype=ComplexF64)  
    mc.thermalization_measurements[:kx] = kx_dia_measurement(mc, model, obs=FullBinner(Float64), eltype=Float64)
    mc.thermalization_measurements[:CDSxx]=charge_density_susceptibility(mc, model, obs=FullBinner(Array{ComplexF64}), 
        flavor_iterator = MonteCarlo.FlavorIterator(mc, 0),   kernel = full_cdc_XX_kernel)

    ##### B1-nematic sector
    mc.thermalization_measurements[:B1_CDS] = B1_charge_density_susceptibility(mc, model, obs=FullBinner(Float64))
    mc.thermalization_measurements[:B1_CDC] = B1_charge_density_correlation(mc, model, obs=FullBinner(Float64))
    mc.thermalization_measurements[:NemC_X] = nematic_correlation(mc, model, kernel =nem_X_kernel, obs=FullBinner(Float64))
    mc.thermalization_measurements[:NemS_X] = nematic_susceptibility(mc, model, kernel =nem_X_kernel, obs=FullBinner(Float64))

    ##### A1p double-Q sector
    mc.thermalization_measurements[:A1p_dQ_C_X] = A1_Q1Q2_correlation(mc, model, kernel=A1_Q1Q2_X_kernel, obs=FullBinner(Float64))
    mc.thermalization_measurements[:A1p_dQ_S_X] = A1_Q1Q2_susceptibility(mc, model, kernel=A1_Q1Q2_X_kernel, obs=FullBinner(Float64))

    ##### B1p double-Q sector
    mc.thermalization_measurements[:B1p_OP_z] = B1p_OP(mc, model, :z, obs=FullBinner(Float64))



    ##############
    # Adding measurements
    ##############
    cap=16*ceil(Int, sweeps /mc.parameters.measure_rate);

    # mc.measurements[:total_greens] = total_greens_measurement(mc, model, capacity=cap)

    ##### energy and heat capacity
    #
    En_eltype=Float64;
    mc.measurements[:E] = total_energy(mc, model, obs=FullBinner(En_eltype))
    mc.measurements[:E_TI] = total_energy(mc, model, greens_iterator = TimeIntegral(mc),
         lattice_iterator = energy_iterator(), obs=FullBinner(En_eltype))

    mc.measurements[:E_kin] = kinetic_energy(mc, model, obs=FullBinner(En_eltype))
    mc.measurements[:E_pot] = interaction_energy(mc, model, obs=FullBinner(En_eltype))
    heat_eltype=Float64;
    
    mc.measurements[:heat_cap_h2] = heat_cap_h2(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.measurements[:heat_cap_h3] = heat_cap_h3(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.measurements[:heat_cap_h4] = heat_cap_h4(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)

    mc.measurements[:heat_cap_h2_TI] = heat_cap_h2_TI(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.measurements[:heat_cap_h3_TI] = heat_cap_h3_TI(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.measurements[:heat_cap_h4_TI] = heat_cap_h4_TI(mc, model, obs=FullBinner(heat_eltype), eltype=heat_eltype)

    mc.measurements[:heat_cap_h4_onsite] = heat_cap_h4(mc, model, lattice_iterator=EachSiteSummed(), obs=FullBinner(heat_eltype), eltype=heat_eltype)
    mc.measurements[:heat_cap_h4_onsite_TI] = heat_cap_h4_TI(mc, model, lattice_iterator=EachSiteSummed(), obs=FullBinner(heat_eltype), eltype=heat_eltype)

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

    mc.measurements[:SDC_A1_Mx_x] = spin_density_correlation(mc, model, :x, kernel=full_sdc_Mx_x_kernel, 
        lattice_iterator = EachSitePair_A1(), capacity=cap)
    mc.measurements[:SDS_A1_Mx_x] = spin_density_susceptibility(mc, model, :x, kernel=full_sdc_Mx_x_kernel, 
        lattice_iterator = EachSitePair_A1(), capacity=cap)



    ##### superconductivity
    #
    mc.measurements[:PDC_s] = pairing_correlation(mc, model, kernel = pc_swave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_s] = pairing_susceptibility(mc, model, kernel = pc_swave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDC_spm] = pairing_correlation(mc, model, kernel = pc_spm_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_spm] = pairing_susceptibility(mc, model, kernel = pc_spm_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_XX] = pairing_susceptibility(mc, model, kernel = pc_XX_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)
    mc.measurements[:PDS_YYzz] = pairing_susceptibility(mc, model, kernel = pc_YYzz_wave_kernel, lattice_iterator = EachSitePairByDistance(), capacity=cap)

    mc.measurements[:Δ_Zy_bil_OP] = Δ_Zy_bil_OP(mc, model)
    mc.measurements[:Δ_0y_bil_OP] = Δ_0y_bil_OP(mc, model)
    mc.measurements[:Δ_Xy_bil_OP] = Δ_Xy_bil_OP(mc, model)
    mc.measurements[:Δ_Ysum_bil_OP] = Δ_Ysum_bil_OP(mc, model)


    ##### charge sector
    #
    #mc.measurements[:occ] = occupation(mc, model, capacity=cap)
    mc.measurements[:occ_summed] =occupation_summed(mc, model, eltype=Float64)
    #mc.measurements[:CDC] = charge_density_correlation(mc, model, capacity=cap)
    mc.measurements[:CDC_summed] = charge_density_correlation(mc, model, lattice_iterator = EachSitePair_summed(), eltype =Float64, obs=FullBinner(Float64))

    mc.measurements[:CDS] = charge_density_susceptibility(mc, model, capacity=cap)
    #mc.measurements[:Zk_proxy] = spectral_weight_proxy(mc, model, capacity=cap)
    # mc.measurements[:green]= greens_measurement(mc, model, capacity=cap)
    # mc.measurements[:CCS]= current_current_susceptibility(mc, model, flavor_iterator = MonteCarlo.FlavorIterator(mc, 0),
    #         kernel = my_cc_kernel, capacity=cap)

    mc.measurements[:rho_s]  = phase_stiffness(mc, model, capacity=cap, eltype=ComplexF64)  
    mc.measurements[:kx] = kx_dia_measurement(mc, model, capacity=cap, eltype=Float64)

    #mc.measurements[:CD_X_OP] = CD_X_OP(mc, model)  #For reasons unclear, the ρ^{x0} -OP is purely imaginary
    #CD_X OP is zero in case of Ising-X
    mc.measurements[:CDSxx]=charge_density_susceptibility(mc, model, capacity=cap, flavor_iterator = MonteCarlo.FlavorIterator(mc, 0),
            kernel = full_cdc_XX_kernel)



    ##### B1-nematic sector
    #
    mc.measurements[:B1_OP] = nematic_OP(mc, model, obs=FullBinner(Float64))
    mc.measurements[:B1_proxy_OP] = proxy_B1_OP(mc, model)
    mc.measurements[:B1_CDS] = B1_charge_density_susceptibility(mc, model, capacity=cap)
    mc.measurements[:B1_CDC] = B1_charge_density_correlation(mc, model, capacity=cap)


    # mc.measurements[:NemS] = nematic_susceptibility(mc, model, capacity=cap)
    #mc.measurements[:NemC] = nematic_correlation(mc, model, capacity=cap)
    mc.measurements[:NemC_X] = nematic_correlation(mc, model, kernel =nem_X_kernel, capacity=cap)
    mc.measurements[:NemS_X] = nematic_susceptibility(mc, model, kernel =nem_X_kernel, capacity=cap)

    ##### A1p double-Q sector
    #
    mc.measurements[:A1p_OP] = A1p_OP(mc, model, obs=FullBinner(Float64))
    mc.measurements[:A1p_proxy_OP] = proxy_A1p_OP(mc, model)
    # mc.measurements[:A1p_dQ_S] = A1_Q1Q2_susceptibility(mc, model, capacity=cap)
    # mc.measurements[:A1p_dQ_C] = A1_Q1Q2_correlation(mc, model, capacity=cap)
    mc.measurements[:A1p_dQ_C_X] = A1_Q1Q2_correlation(mc, model, kernel=A1_Q1Q2_X_kernel, capacity=cap)
    mc.measurements[:A1p_dQ_S_X] = A1_Q1Q2_susceptibility(mc, model, kernel=A1_Q1Q2_X_kernel, capacity=cap)

    ##### B1p double-Q sector
    #
    #mc.measurements[:B1p_OP_x] = B1p_OP(mc, model, :x, obs=FullBinner(Float64))
    #mc.measurements[:B1p_OP_y] = B1p_OP(mc, model, :y, obs=FullBinner(Float64))
    mc.measurements[:B1p_OP_z] = B1p_OP(mc, model, :z, obs=FullBinner(Float64))
    #B1p OP is zero in x, and y direction in case of Ising-X
    # mc.measurements[:B1p_dQ_C] = B1p_Q1Q2_correlation(mc, model, capacity=cap)
    # mc.measurements[:B1p_dQ_S] = B1p_Q1Q2_susceptibility(mc, model, capacity=cap)


    ################
    # run the simulation
    ###############


    exit_code= run!(mc, sweeps=sweeps, thermalization=therm, verbose=true,  ten_reg=3, _measure_th_config=true)

    # path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/projects/IsingX_Striped_d/materials/"
    # st="IsingX_d_save_load_test.jld2";
    # final_file=path * "/" * st;

    # MonteCarlo.save(final_file, mc, overwrite = true, rename = false)

