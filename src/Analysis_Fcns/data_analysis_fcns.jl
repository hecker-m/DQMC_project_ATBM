using Interpolations, Roots
γOν=1.5;
μ0=1.5;
T_max=0.4;
interp_linear=Vector(undef, 4)
crossing_fcns=Vector{Function}(undef, factorial(3))
corssing_values=Vector{Float64}(undef, factorial(3))
for (idx, Ls) in enumerate(unique(df_LT[:,:L]))
    dfs=filter(:T => T -> T≤T_max  ,filter(:L=>L -> L==Ls ,filter(:μ0 =>μ ->μ==μ0, dfU0)))  

    global listx=dfs[!,:T]
    global listy=dfs[!,Symbol("SDS_A1_Mx_x")]  .*(Ls^(2-γOν))
    interp_linear[idx] = linear_interpolation(dfs[!,:T], dfs[!,Symbol("SDS_A1_Mx_x")]  .*(Ls^(2-γOν)))

end
tot=1;
for i in range(1,4)
    for j in range(i+1,4)
        global tot
        func(x)=interp_linear[i](x)-interp_linear[j](x)
        crossing_fcns[tot]=func
        corssing_values[tot]=find_zero(func, (0.1, 0.4))
        tot+=1
    end
end

"""
    compute_crossing(x_list, y_lists)
return the approximate x-value at which the curves in `y_lists` intersect 
"""
function compute_crossing(x_list, y_lists)
    N_y=length(y_lists)
    xmin=minimum(x_list)
    xmax=maximum(x_list)

    interp_linear=Vector{Any}(undef, N_y)
    crossing_values=Vector{Float64}(undef, factorial(N_y-1))
    for idx in 1:N_y
        interp_linear[idx] = linear_interpolation(x_list, y_lists[idx])
    end
    tot=1;
    for i in range(1,N_y)
        for j in range(i+1,N_y)
            func(x)=interp_linear[i](x)-interp_linear[j](x)
            crossing_values[tot]=find_zero(func, (xmin, xmax))
            tot+=1
        end
    end
    return mean(crossing_values)
end


function compute_crossing_γOν_full(γOν, x_list, χ_lists, L_list)
    return compute_crossing(x_list, χ_lists .*(L_list .^(-γOν)))
end


Spin_scaling_γOν_data=DataFrame(U=Float64[], μ0=Float64[], γOν=Float64[])


U_list=[1.75, ]
U_list=[1.25, 1.5, 1.75, 2.0, 2.25]

for U0 in U_list
    dfU0=df_U[df_U.keymap[(U0,)]]
    df_LT=sort(unique(dfU0[:,[:T,:L]]),:T, rev=true)
    for μ0 in unique(df_Lμ[:,:μ0])
    # for μ0 in [1.3, ]

        key=(U0, μ0);
        if haskey(T_spin_grouped.keymap, key)
            global x_list=filter(:T => T -> T≤T_max  ,filter(:L=>L -> L==unique(df_LT[:,:L])[1] ,
                filter(:μ0 =>μ ->μ==μ0, dfU0)))[!,:T]
            ;
            global χ_lists=[filter(:T => T -> T≤T_max  ,filter(:L=>L -> L==Ls ,
                filter(:μ0 =>μ ->μ==μ0, dfU0)))[!,Symbol("SDS_A1_Mx_x")]  .*(Ls^(2)) 
                for Ls in unique(df_LT[:,:L])]
            ;    
            L_list=[Ls for Ls in unique(df_LT[:,:L])]
            # @show L_list
            # @show χ_lists
            Tc_spin=T_spin_grouped[T_spin_grouped.keymap[key]].Tc_spin[1]

            global compute_crossing_γOν(γOν)=compute_crossing_γOν_full(γOν, x_list, χ_lists, L_list)-Tc_spin


            function f_γOν_max(γOν)
                return χ_lists[1][1]*(L_list[1]^(-γOν)) -χ_lists[4][1]*(L_list[4]^(-γOν))
            end
            # vec=[f_γOν_max(γOν) for γOν in range(1,2,5)]

            γOν_max=find_zero(f_γOν_max,1.5)

            function f_γOν_min(γOν)
                NT=length(χ_lists[1])
                return χ_lists[1][NT]*(L_list[1]^(-γOν)) -χ_lists[4][NT]*(L_list[4]^(-γOν))
            end
            # vec=[f_γOν_min(γOν) for γOν in range(0.1,1,5)]

            γOν_min=find_zero(f_γOν_min, 0.1)

            @show (γOν_min, γOν_max)

            # γOν_value=find_zero(compute_crossing_γOν, (γOν_min, γOν_max))
            local γOν_value
            try
                γOν_value=fzero(compute_crossing_γOν, γOν_min+(γOν_max-γOν_min)*9/10)
            catch
                println("!!! U=$(U0), μ=$(μ0), has not worked. γOν_min=$(γOν_min), γOν_max=$(γOν_max)")

            else
                push!(Spin_scaling_γOν_data, [U0, μ0, γOν_value])
                println("done with U=$(U0), μ=$(μ0)")
            end

        end
    end
end



#####################
### Plotting the derived γ/ν values
####################
Spin_scaling_γOν_data_U=groupby(Spin_scaling_γOν_data, :U, sort=true)


fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$μ/t$", ylabel=L"$γ/ν$", ylabelsize=30,
    xlabelsize=30
)

for (idx, U0) in enumerate(unique(Spin_scaling_γOν_data[:,:U]))
    γOν_data=Spin_scaling_γOν_data_U[Spin_scaling_γOν_data_U.keymap[(U0,)]]

    CairoMakie.scatterlines!(top, γOν_data[:, :μ0], γOν_data[:, :γOν] ;  
        marker =:xcross, markersize=15,  linewidth=2,
        color = colorschemes[:tab20][2*idx-1], label = L"$U=%$(U0)$"
    )
    
end
axislegend( position=(1, 1))

display(fig)
if save_bool_fig
    global pPap=p * "paper_pics/" 
    !isdir(pPap) ? mkdir(pPap) : nothing ;
    CairoMakie.save(joinpath(pPap, "Gamma_over_Nu"  * ".pdf"), fig; 
        pt_per_unit=2 
    )
end


