################
## loading functions
################
"""
Similar to the MonteCarlo.load() function, but 

my_load(data, ::Val{:DQMC}; _recorder=true, _th_meas=true, _meas=true)

saves timebut optionally not loading e.g. thermalization measurements, 
or the recorded configurations
"""
function my_load(data, ::Val{:DQMC}; _recorder=true, _th_meas=true, _meas=true)
    if data["VERSION"] > 3
        throw(ErrorException("Failed to load DQMC version $(data["VERSION"])"))
    end
    parameters = MonteCarlo._load(data["Parameters"])
    if haskey(data, "CB")
        parameters = DQMCParameters(parameters, checkerboard = data["CB"])
    end
    tag = Val(Symbol(get(data["Analysis"], "tag", :DQMCAnalysis)))
    analysis = MonteCarlo._load(data["Analysis"], tag)
    last_sweep = data["last_sweep"]
    model = MonteCarlo.load_model(data["Model"], MonteCarlo.to_tag(data["Model"]))

    if _recorder
        recorder = MonteCarlo._load(data["configs"], MonteCarlo.to_tag(data["configs"]))
    else
        #field = MonteCarlo.choose_field(model);
        recorder = ConfigRecorder(AbstractMagnBosonField, 10) 
        #initialize empty recorder to save memory
    end

    if haskey(data, "field")
        tag = Val(Symbol(get(data["field"], "tag", :Field)))
        field = MonteCarlo._load(data["field"], tag, parameters, model)
    else
        conf = data["conf"]
        field = field_hint(model, to_tag(data["Model"]))(parameters, model)
        conf!(field, conf)
    end
    scheduler = if haskey(data, "Scheduler")
        MonteCarlo._load(data["Scheduler"])
    else
        if haskey(data["Parameters"], "global_moves") && Bool(data["Parameters"]["global_moves"])
            rate = get(data["Parameters"], "global_rate", 10)
            @warn "Replacing `global_moves = true` with GlobalFlip"
            SimpleScheduler(LocalSweep(rate), GlobalFlip())
        else
            SimpleScheduler(LocalSweep())
        end
    end


   #tag = Val(Symbol(data["Measurements"], "tag", :Measurements))
   #I do not know what tag was good for, but it is not defined well!!!
   #Loads the whole thing, and slows the shit down.
    if _th_meas && _meas
        combined_measurements = MonteCarlo._load(data["Measurements"])
        thermalization_measurements = combined_measurements[:TH]
        measurements = combined_measurements[:ME]
    elseif !_th_meas && _meas
        thermalization_measurements = Dict{Symbol, AbstractMeasurement}()
        measurements = MonteCarlo._load(data["Measurements"])[:ME]
    elseif _th_meas && !_meas
        thermalization_measurements = MonteCarlo._load(data["Measurements"])[:TH]
        measurements = Dict{Symbol, AbstractMeasurement}()
    else
        thermalization_measurements = Dict{Symbol, AbstractMeasurement}()
        measurements = Dict{Symbol, AbstractMeasurement}()
    end

    stack = DQMCStack(field, model, parameters.checkerboard)
    ut_stack = MonteCarlo.UnequalTimeStack{MonteCarlo.geltype(stack), MonteCarlo.gmattype(stack)}()

    DQMC(
        model, field, last_sweep, 
        stack, ut_stack, scheduler,
        parameters, analysis, 
        recorder, thermalization_measurements, measurements
    )
end

"""
To overcome conflicts in MonteCarlo.jl versioning,
I defined my own load function

_load(dqmcs::Array, L::Int, T::Real, beta::Real, U::Real, peierls::Bool, 
    therm::Int, sweeps::Int, N_worker::Int, beta_bool::Bool, jobid::Int=0; 
   path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/Reprod_Fernandes_Paper/run_saves/", 
   prefix_folder="", prefix_file="FP_", kwargs...)

As a consistency check it returns the number of workers whose DQMC simulations 
were loaded into the array dqmcs.

kwargs are _recorder=true, _th_meas=true, _meas=true from my_load(data["MC"], Val(:DQMC); kwargs...)
"""
function _load(dqmcs::Array, L::Int, T::Real, beta::Real, U::Real, peierls::Bool, 
     therm::Int, sweeps::Int, N_worker::Int, beta_bool::Bool, jobid::Int=0; 
    path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/Reprod_Fernandes_Paper/run_saves/", 
    prefix_folder="", prefix_file="FP_", eps="", kwargs...)
    
    beta_bool ? tempstring="b_" * to_string(beta) : tempstring="T_" * to_string(T)

    n_workers=0;
    for worker =1:N_worker
        _path=path * prefix_folder * "L$(L)_" * tempstring * "_U_" * to_string(U) * eps * "_B_$(Int(peierls))";
        jobid==0 && begin _filename=prefix_file * tempstring *"_U_"* to_string(U) *
            "_L$(L)" * eps * "_B_$(Int(peierls))_sw$(sweeps)_th$(therm)_worker_$(worker).jld2";
        end 
        jobid==0 || begin _filename=prefix_file * tempstring *"_U_"* to_string(U) *
        "_L$(L)" * eps * "_B_$(Int(peierls))_sw$(sweeps)_th$(therm)_worker_$(worker)_id$(jobid).jld2";
        end 

        filename=_path * "/" * _filename;
        if isfile(filename)
            data = MonteCarlo.FileData(JLD2.load(filename), abspath(filename));
            mc=my_load(data["MC"], Val(:DQMC); kwargs...)
            push!(dqmcs, mc);
            n_workers+=1;
        end
    end  
    return n_workers
end
#
#_load_full can be used. I only need to comment 
#tag = Val(Symbol(data["Measurements"], "tag", :Measurements))
#(see above), inside the load function.
"""
_load_full(dqmcs, L::Int, T::Float64, U::Float64, peierls::Bool, 
    jobid::Int, therm::Int, sweeps::Int, N_worker::Int; prefix::String="FP_")

loads the full DQMC simulation, found at the specified path. 

"""
function _load_full(dqmcs, L::Int, T::Float64, U::Float64, peierls::Bool, 
    jobid::Int, therm::Int, sweeps::Int, N_worker::Int; prefix::String="FP_", eps="")
    for worker =1:N_worker
        path="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/" *
            "Reprod_Fernandes_Paper/run_saves/L$(L)_T_" *
            to_string(T) * "_U_" * to_string(U) * eps * "_B_$(Int(peierls))";
        st=prefix*"T_" * to_string(T) *"_U_"* to_string(U) *
        "_L$(L)" * eps * "_B_$(Int(peierls))_sw$(sweeps)_th$(therm)_worker_$(worker)_id$(jobid).jld2";
        filename=path * "/" * st;
        mc=MonteCarlo.load(filename);
        push!(dqmcs, mc);
    end  
end
