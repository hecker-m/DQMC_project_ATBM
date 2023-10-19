function binder_fcn(S_i, S2_i)
    return  1 - S2_i / (3*S_i^2)
end


function calculate_Binder(S_s::Vector , S2_s::Vector , bin_length::Integer)

    binned_S, binned_S2 = prebinning(bin_length, S_s , S2_s)

    #μ_g, σJackKnife_g = calc_mean_JackKnifeError(binder_fcn, binned_S, binned_S2)
    μ_g, σJackKnife_g = jackknife(binder_fcn, binned_S, binned_S2)
    return μ_g, σJackKnife_g
end


"""
prebinning(bin_length::Integer, vecs::Vector{<:Number}...)
collapses the long vectors in `vecs` into binned vectors 
where each bin comprises `bin_length` number of former values.
"""
function prebinning(bin_length::Integer, vecs::Vector{<:Number}...)
    if mod(length(vecs[1]), bin_length) != 0
        println("Assumed binner length $(bin_length) is no integer multiple of data length $(length(vecs[1])).
                Left $(mod(length(vecs[1]), bin_length)) points unaccounted.")
    end
    N_bins=div(length(vecs[1]), bin_length)

    binned_vectors = map(vecs) do vec
        [mean(vec[1 + (n-1)*bin_length : n*bin_length] ) for n in 1:N_bins]
    end

    return binned_vectors
end




"""
    calc_mean_JackKnifeError(g::Function, x::Vector, y::Vector)
    computes the JackKnife estimate and the corresponding error.
    Has become obsolete since BinningAnalysis.jackknife() does the same,
    and for arbitrary number of arguments and complex values.  
"""
function calc_mean_JackKnifeError(g::Function, x::Vector, y::Vector)
    M=length(x)
    x_bar =mean(x)
    y_bar = mean(y)
    g_xy_bar= g(x_bar, y_bar)
    g_i=Vector{Float64}(undef, M)
    for i in 1:M
        g_i[i]= g(x_bar + (x_bar -x[i])/(M-1),  y_bar + (y_bar -y[i])/(M-1))
    end
    g_estimator = M*g_xy_bar -(M-1) * mean(g_i)
    g_JK_error = sqrt((M-1)/M * sum( (g_i[:] .- g_xy_bar) .^2 ))

    return g_estimator, g_JK_error
end