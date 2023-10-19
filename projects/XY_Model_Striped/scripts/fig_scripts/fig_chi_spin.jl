

########################
## spin susceptibility, Q=(0,0)
########################
q0=["0", "0"]
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}=(0,0)]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=$(L_plot)")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("SDS_Mx_x" * "_$(q0[1])$(q0[2])")] .+dfs[!,Symbol("SDS_Mx_y" * "_$(q0[1])$(q0[2])")], 
            dfs[!,Symbol("ΔSDS_Mx_x" * "_$(q0[1])$(q0[2])")] .+dfs[!,Symbol("ΔSDS_Mx_y" * "_$(q0[1])$(q0[2])")]; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDS_Mx_x" * "_$(q0[1])$(q0[2])")] .+dfs[!,Symbol("SDS_Mx_y" * "_$(q0[1])$(q0[2])")] ;  marker = :xcross,
            color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=20)


    CairoMakie.lines!(top, range(Umin, Umax,nU), [2*χ_spin_XX(1/Ts, U/Nϕ, muM) for U in range(Umin, Umax,nU)] ;  
            linestyle = :dash,   color = colorschemes[:coolwarm][1-1/Ts/βmax])

end
axislegend( position=(0.04,0.97))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q($(q0[1]),$(q0[2])).png"), fig)

########################
## spin susceptibility, Q₁ =(π,0) + Q₂=(0,π)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}_1] +
    χ_{\mathrm{spin}}^{} [\mathbf{Q}_2]$", ylabelsize=30,
    xlabelsize=30, title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")] .+
            dfs[!,Symbol("SDS_Mx_y_0π")] .+ dfs[!,Symbol("SDS_Mx_y_π0")], 
            dfs[!,Symbol("ΔSDS_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_x_π0")] .+dfs[!,Symbol("ΔSDS_Mx_y_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_y_π0")],
            df_ϕs[!,Symbol("ΔSDS_A1_00")]*Ls^2; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")].+dfs[!,Symbol("SDS_Mx_y_0π")] .+ dfs[!,Symbol("SDS_Mx_y_π0")] ;
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("SDS_A1_00")]*Ls^2 , 
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax], 
         markersize=15)     
         
    CairoMakie.lines!(top, range(Umin, Umax,nU), [4*χ_spin_XX(1/Ts, U/Nϕ, muM) for U in range(Umin, Umax,nU)] ;  
         linestyle = :dash,   color = colorschemes[:coolwarm][1-1/Ts/βmax])
end

axislegend( position=(0, 1))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q_0π_Plus_π0.png"), fig)

########################
## spin susceptibility, Q₁=(π,0) + Q₂=(0,π), as function of temperature
########################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1], xticklabelsvisible=false, ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}_1] \!+\!\!
    χ_{\mathrm{spin}}^{} [\mathbf{Q}_2]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=$(L_plot)")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  


    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")] .+dfs[!,Symbol("SDS_Mx_y_0π")] .+ dfs[!,Symbol("SDS_Mx_y_π0")], 
            dfs[!,Symbol("ΔSDS_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_x_π0")] .+dfs[!,Symbol("ΔSDS_Mx_y_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_y_π0")],
            df_ϕs[!,Symbol("ΔSDS_A1_00")]*Ls^2; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")] .+dfs[!,Symbol("SDS_Mx_y_0π")] .+ dfs[!,Symbol("SDS_Mx_y_π0")] ;
        marker = :xcross, color = colorschemes[:tab20][1-idx/num_T_points], 
        label = "U=$(Us)", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("SDS_A1_00")]*Ls^2 , 
        marker = '□', color = colorschemes[:tab20][1-idx/num_T_points], 
         markersize=15)         
end
axislegend( position=(1, 1))

########################
## ZZ - spin susceptibility, Q₁=(π,0) + Q₂=(0,π), as function of temperature
########################
top=Axis(fig[2, 1], xlabel=L"$T/t$", ylabel=L"$χ_{\mathrm{spin}}^{z} [\mathbf{Q}_1] \!+\!\!
    χ_{\mathrm{spin}}^{z} [\mathbf{Q}_2]$", ylabelsize=30,
    xlabelsize=30)

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  


    eb=errorbars!(dfs[!,:T], dfs[!,Symbol("SDS_Mx_z_0π")] .+ dfs[!,Symbol("SDS_Mx_z_π0")], 
            dfs[!,Symbol("ΔSDS_Mx_z_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_z_π0")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], dfs[!,Symbol("SDS_Mx_z_0π")] .+ dfs[!,Symbol("SDS_Mx_z_π0")]  ;
        marker = :xcross, color = colorschemes[:tab20][1-idx/num_T_points], 
        label = "U=$(Us)", markersize=10)
       
end

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q_0π_Plus_π0_T.png"), fig)

################
## spin correlation function, Q₁=(π,0) + Q₂=(0,π)
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$\mathrm{corr:\,}S_{\mathrm{spin}} [\mathbf{Q}_1] \!+\!\!
    S_{\mathrm{spin}} [\mathbf{Q}_2]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ))  

    eb=errorbars!(dfs[!,:U], (dfs[!,Symbol("SDC_Mx_x_0π")] .+ dfs[!,Symbol("SDC_Mx_x_π0")] 
        .+dfs[!,Symbol("SDC_Mx_y_0π")] .+ dfs[!,Symbol("SDC_Mx_y_π0")]) /Ls^2, 
        (dfs[!,Symbol("ΔSDC_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDC_Mx_x_π0")] .+dfs[!,Symbol("ΔSDC_Mx_y_0π")] .+ dfs[!,Symbol("ΔSDC_Mx_y_π0")]) /Ls^2, 
        df_ϕs[!,Symbol("ΔSDC_A1_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], (dfs[!,Symbol("SDC_Mx_x_0π")] .+ dfs[!,Symbol("SDC_Mx_x_π0")] .+dfs[!,Symbol("SDC_Mx_y_0π")] .+ dfs[!,Symbol("SDC_Mx_y_π0")]) /Ls^2  ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:U], df_ϕs[!,Symbol("SDC_A1_00")] ;  
            marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15) 

    CairoMakie.lines!(top, range(Umin, Umax,nU), [4/Ls^2*S_spin_XX(1/Ts, U/Nϕ, muM) for U in range(Umin, Umax,nU)] ;  
            linestyle = :dash,   color = colorschemes[:coolwarm][1-1/Ts/βmax])
end
#axislegend( position=(0, 1))


########################
## Magnetic Order parameter
########################
top=Axis(fig[2, 1], xlabel=L"$U/t$", ylabel=L"$\langle \left| \left( \mathbf{M}_{\mathbf{Q}_1}, 
    \mathbf{M}_{\mathbf{Q}_2}   \right)  \right| \rangle $", ylabelsize=30,
    xlabelsize=30, limits= (nothing, nothing))

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), (0, 0, -0.8))

for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP)) 
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_OP))  
 
    eb=errorbars!(df_OPs[!,:U], df_OPs[!,Symbol("Mx_X_OP34")], df_OPs[!,Symbol("ΔMx_X_OP34")],
        df_ϕ_OPs[!,Symbol("ΔMx_OP")]/sqrt(2) ; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:U], df_OPs[!,Symbol("Mx_X_OP34")] ;  marker = :xcross,
        color = colorschemes[:coolwarm][1-1/Ts/βmax], label = "β=$((1/Ts))", markersize=10) 
        
    CairoMakie.scatter!(top, df_ϕ_OPs[!,:U], df_ϕ_OPs[!,Symbol("Mx_OP")]/sqrt(2)  ;  
        marker = '□', color = colorschemes[:coolwarm][1-1/Ts/βmax],  markersize=15)         
end

display(fig)
CairoMakie.save(joinpath(p, "corr_spin.png"), fig)



################
## spin correlation function, Q=(π,0) + Q=(0,π), as a function of T
################
fig = Figure(resolution = (800, 800))
top=Axis(fig[1, 1],  ylabel=L"$ S_{\mathrm{spin}} [\mathbf{Q}_1] \!+\!\! 
    S_{\mathrm{spin}} [\mathbf{Q}_2]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing),
    title=L"B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(dfs[!,:T], (dfs[!,Symbol("SDC_Mx_x_0π")] .+ dfs[!,Symbol("SDC_Mx_x_π0")] .+dfs[!,Symbol("SDC_Mx_y_0π")] .+ dfs[!,Symbol("SDC_Mx_y_π0")])/Ls^2, 
        (dfs[!,Symbol("ΔSDC_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDC_Mx_x_π0")] .+dfs[!,Symbol("ΔSDC_Mx_y_0π")] .+ dfs[!,Symbol("ΔSDC_Mx_y_π0")]) /Ls^2, 
        df_ϕs[!,Symbol("ΔSDC_A1_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], (dfs[!,Symbol("SDC_Mx_x_0π")] .+ dfs[!,Symbol("SDC_Mx_x_π0")] .+dfs[!,Symbol("SDC_Mx_y_0π")] .+ dfs[!,Symbol("SDC_Mx_y_π0")]) /Ls^2  ;  marker = :xcross,
        color = color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:T], df_ϕs[!,Symbol("SDC_A1_00")] ;  
            marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=15) 
end
axislegend( position=(1, 1))
################
## ZZ - spin correlation function, Q=(π,0) + Q=(0,π), as a function of T
################
top=Axis(fig[2, 1],  ylabel=L"$ S^z_{\mathrm{spin}} [\mathbf{Q}_1] \!+\!\! 
    S^z_{\mathrm{spin}} [\mathbf{Q}_2]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Us, Ls)) in enumerate(zip(df_LU[:,:U],df_LU[:,:L]))
    num_T_points=length(df_LU[:,:U]);
    dfs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ))  

    eb=errorbars!(dfs[!,:T], (dfs[!,Symbol("SDC_Mx_z_0π")] .+ dfs[!,Symbol("SDC_Mx_z_π0")])/Ls^2, 
        (dfs[!,Symbol("ΔSDC_Mx_z_0π")] .+ dfs[!,Symbol("ΔSDC_Mx_z_π0")]) /Ls^2  ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:T], (dfs[!,Symbol("SDC_Mx_z_0π")] .+ dfs[!,Symbol("SDC_Mx_z_π0")]) /Ls^2  ;  marker = :xcross,
        color = color = colorschemes[:tab20][1-idx/num_T_points], label = "U=$(Us)", markersize=10)

end

########################
## Magnetic Order parameter, as a function of temperature
########################
top=Axis(fig[3, 1], xlabel="T/t", ylabel=L"$\left| \left(  \mathbf{M}_{\mathbf{Q}_1}, \!
\mathbf{M}_{\mathbf{Q}_2}  \right)  \right|$", xlabelsize=30, ylabelsize=30,
xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

mkey="Mx_X_OP34"

df_OP_LU=sort(unique(df_OP[:,[:U,:L]]),:U)
for (idx, (Us, Ls)) in enumerate(zip(df_OP_LU[:,:U],df_OP_LU[:,:L]))
    num_T_points=length(df_OP_LU[:,:U]);
    df_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_OP))  
    df_ϕ_OPs=filter(:L=>L -> L==Ls ,filter(:U =>U ->U==Us, df_ϕ_OP))  


    eb=errorbars!(df_OPs[!,:T], df_OPs[!,Symbol(mkey)], df_OPs[!,Symbol("Δ" * mkey)]; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_T_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:T], df_OPs[!,Symbol(mkey)] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_T_points], 
        label = "U=$(round(Us, digits=2)), L=$Ls", markersize=10)  
        
    # CairoMakie.scatter!(top, df_ϕ_OPs[!,:T], df_ϕ_OPs[!,Symbol("Mx_OP")]/2  ;  
    #     marker = '□', color = colorschemes[:tab20][1-idx/num_T_points],  markersize=15)      
end
#axislegend( position=(1,1))

display(fig)
CairoMakie.save(joinpath(p, "corr_spin_T.png"), fig)




########################
## spin susceptibility, Q₁=(π,0) - Q₂=(0,π)
########################
fig = Figure(resolution = (800, 600))
top=Axis(fig[1, 1], xlabel=L"$U/t$", ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}_1] \!-\!\!
    χ_{\mathrm{spin}}^{} [\mathbf{Q}_2]$", ylabelsize=30,
    xlabelsize=30, title="B=$(Int(peierls)),  L=$(L_plot)")

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_U_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df))  
    eb=errorbars!(dfs[!,:U], dfs[!,Symbol("SDS_Mx_x_0π")] .- dfs[!,Symbol("SDS_Mx_x_π0")] .+dfs[!,Symbol("SDS_Mx_y_0π")] .- dfs[!,Symbol("SDS_Mx_y_π0")], 
            dfs[!,Symbol("ΔSDS_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_x_π0")] .+dfs[!,Symbol("ΔSDS_Mx_y_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_y_π0")]; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:coolwarm][1-1/Ts/βmax])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:U], dfs[!,Symbol("SDS_Mx_x_0π")] .- dfs[!,Symbol("SDS_Mx_x_π0")] .+dfs[!,Symbol("SDS_Mx_y_0π")] .- dfs[!,Symbol("SDS_Mx_y_π0")] ;
        marker = :xcross, color = colorschemes[:coolwarm][1-1/Ts/βmax], 
        label = "β=$((1/Ts))", markersize=10)
end
axislegend( position=(0.04,0.97))

display(fig)
CairoMakie.save(joinpath(p, "chi_spin_Q_0π_Minus_π0.png"), fig)


