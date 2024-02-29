include("/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/Analysis_Fcns/analysis_fcns.jl")
include("/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/Analysis_Fcns/analysis_fcns_phi.jl")




########################################################
########################################################
## Analysis function for thermalization measurements
## Note that thermalization measurement are always done in FullBinner
########################################################
########################################################

function convert_observable(mc, tuple::NamedTuple, _th_measures::Int)
    f_key =first_key(tuple.key)
    return real(mc.thermalization_measurements[f_key].observable.x[_th_measures])
end


function convert_observable_ϕ(mc, _th_measures::Int, ::Val{:Q0πQπ0_offset})
    U=mc.model.U
    lat=mc.model.l
    L=lat.Ls[1]
    N=length(lat)
    Nϕ=mc.parameters.Nϕ
    Nτ=mc.parameters.slices
    δτ=mc.parameters.delta_tau
    β=mc.parameters.beta

    ϕ_field=extract_ϕ(mc, _th_measures)

    ϕQ1Q2= calc_ϕQ1Q2(mc, ϕ_field)
    ΛA1=calc_ΛA1(mc, ϕ_field)
    Ω4A1=calc_ΩnA1(mc, ϕ_field,4)

    ϕQ1Q2_OP=calc_ϕQ1Q2_OP(mc, ϕQ1Q2)
    ϕQ1Q2_bar=calc_ϕQ1Q2_bar(mc, ϕQ1Q2)

    ΦA1=calc_ΦA1(mc, ϕQ1Q2)
    ΦB1=calc_ΦB1(mc, ϕQ1Q2)
    ΦA1p=calc_ΦA1p(mc, ϕQ1Q2)
    ΦB1p=calc_ΦB1p(mc, ϕQ1Q2)

    Epot=-U * (mean(ΛA1) - 1/(2*δτ*U));
    h4=U^2 * (mean(ΛA1 .^2) - (Nϕ+2/N)/(δτ*U) * mean(ΛA1) + Nϕ*(Nϕ+2/N)/(4*δτ^2 *U^2));
    h4_OS=U^2/N * (mean(Ω4A1) - 3/(δτ*U) * mean(ΛA1) +  3/(4*δτ^2 *U^2));

    S_spin=mean(ΦA1) -Nϕ/(N*U*δτ);
    S2_spin=mean(ΦA1 .^2)  - 2*(1+Nϕ)/(N*U*δτ) * mean(ΦA1) + Nϕ*(1+Nϕ)/(δτ^2*N^2 *U^2);
    χ_spin=1/(2U) *(β* sum(ϕQ1Q2_bar .^2) - 2Nϕ /(N));

    S_bil_B1=mean(ΦB1 .^2) - 2/(N*U*δτ) * mean(ΦA1) + Nϕ/(δτ^2*N^2 *U^2);
    S2_bil_B1=mean(ΦB1 .^4)  -12/(N* δτ *U ) *mean(ΦB1 .* ΦB1 .* ΦA1) +
        6*(Nϕ +4 )/(N^2* δτ^2 *U^2 ) *mean(ΦB1 .* ΦB1) + 12/(N^2* δτ^2 *U^2 ) *mean(ΦA1 .* ΦA1) -
        12*(Nϕ +2 )/(N^3* δτ^3 *U^3 ) *mean(ΦA1) + 3Nϕ*(Nϕ +2 )/(N^4* δτ^4 *U^4 );
    χ_bil_B1=β*mean(ΦB1)^2 - 2/(N*U) * mean(ΦA1) + Nϕ/(δτ*N^2 *U^2);

    S_bil_A1p=mean(ΦA1p .^2) - 2/(N*U*δτ) * mean(ΦA1) + Nϕ/(δτ^2*N^2 *U^2);
    S2_bil_A1p=mean(ΦA1p .^4)  -12/(N* δτ *U ) *mean(ΦA1p .* ΦA1p .* ΦA1) +
        6*(Nϕ +4 )/(N^2* δτ^2 *U^2 ) *mean(ΦA1p .* ΦA1p) + 12/(N^2* δτ^2 *U^2 ) *mean(ΦA1 .* ΦA1) -
        12*(Nϕ +2 )/(N^3* δτ^3 *U^3 ) *mean(ΦA1) + 3Nϕ*(Nϕ +2 )/(N^4* δτ^4 *U^4 );
    χ_bil_A1p=β*mean(ΦA1p)^2 - 2/(N*U) * mean(ΦA1) + Nϕ/(δτ*N^2 *U^2);

    return Epot*N, h4, h4_OS, S_spin, S2_spin, χ_spin, S_bil_B1, S2_bil_B1, 
        χ_bil_B1, S_bil_A1p, S2_bil_A1p, χ_bil_A1p
end
