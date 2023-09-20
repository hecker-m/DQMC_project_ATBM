function binder_fcn(S_i, S2_i)
    return  1 - S2_i / (3*S_i^2)
end

function calculate_Binder(S_s::Vector , S2_s::Vector , bin_length::Integer)

    binned_S, binned_S2 = prebinning(S_s , S2_s , bin_length)

    #μ_g, σJackKnife_g = calc_mean_JackKnifeError(binder_fcn, binned_S, binned_S2)
    μ_g, σJackKnife_g = jackknife(binder_fcn, binned_S, binned_S2)

end

"""
prebinning(x_vec::Vector , y_vec::Vector , bin_length::Integer)
collapses the long vectors x,y into binned vectors 
where each bin comprises `bin_length` number of former values.
"""
function prebinning(x_vec::Vector , y_vec::Vector , bin_length::Integer)
    if mod(length(x_vec), bin_length) != 0
        println("Assumed binner length $(bin_length) is no integer multiple of data length $(length(x_vec)).
                Left $(mod(length(x_vec), bin_length)) points unaccounted.")
    end
    N_bins=div(length(x_vec), bin_length)

    binned_x=Vector{Float64}(undef, N_bins)
    binned_y=Vector{Float64}(undef, N_bins)
    for n in 1:N_bins
        binned_x[n] = mean(x_vec[1 + (n-1)*bin_length : n*bin_length] )
        binned_y[n] = mean(y_vec[1 + (n-1)*bin_length : n*bin_length] )
    end

    return binned_x, binned_y
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