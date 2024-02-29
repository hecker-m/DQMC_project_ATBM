include("../../../../src/MonteCarlo.jl_modified/src/my_files/helpful_fcns.jl")


df_UTL=sort(unique(df[:,[:U,:T,:L]]),[:U,:T], rev=[false, true])
df_OP_UTL=sort(unique(df_OP[:,[:U,:T,:L]]),[:U,:T], rev=[false, true])
df_ϕ_UTL=sort(unique(df_ϕ[:,[:U,:T,:L]]),[:U,:T], rev=[false, true])
df_ϕ_OP_UTL=sort(unique(df_ϕ_OP[:,[:U,:T,:L]]),[:U,:T], rev=[false, true])




for U0 in [0.8, 0.9, 1.0, 1.1]
################
## B1 nematic susceptibility as a function of strain ϵ₀
################
fig = Figure(resolution = (800, 1000))
top=Axis(fig[1, 1],  ylabel=L"$χ_{\mathrm{nematic}}^{B_1} [\mathbf{Q}_0]$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= ((-0.005, 0.25), nothing),
    title=L"U=%$(U0),\; B=%$(Int(peierls)),\;  L=%$(L_plot):\quad \times \;↔\;  
    \langle \hat{c}\hat{c}.. \rangle, \qquad □ \;↔\;  \langle ϕ.. \rangle")

CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.02], color = :gray,linewidth=0.2), 0, 0, -0.8)


df_U0=filter(:U =>U ->U==U0, df)
df_ϕ_U0=filter(:U =>U ->U==U0, df_ϕ)
df_OP_U0=filter(:U =>U ->U==U0, df_OP)
df_ϕ_OP_U0=filter(:U =>U ->U==U0, df_ϕ_OP)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T], df_LT[:,:L]))
    num_points=length(unique(df_LT[:,:T]));
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_U0)) 
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_U0)) 

    eb=errorbars!(dfs[!,:ϵ0], dfs[!,Symbol("NemS_00")] , 
            dfs[!,Symbol("ΔNemS_00")], df_ϕs[!,Symbol("ΔNemS_00")] ; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:ϵ0], dfs[!,Symbol("NemS_00")]  ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_points], label = "β=$((1/Ts))", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:ϵ0], df_ϕs[!,Symbol("NemS_00")];  
        marker = '□', color = colorschemes[:tab20][1-idx/num_points],  markersize=15) 
               
end
axislegend( position=(1, 1))

################
## Magnetic susceptibility as a function of strain ϵ₀
################
top=Axis(fig[2, 1], xticklabelsvisible=false, ylabel=L"$χ_{\mathrm{spin}}^{} [\mathbf{Q}_1] \! +\!\!
    χ_{\mathrm{spin}}^{} [\mathbf{Q}_2]$", ylabelsize=30, limits= ((-0.005, 0.25), nothing),
    xlabelsize=30)

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_U0))  
    df_ϕs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_ϕ_U0))  

    eb=errorbars!(dfs[!,:ϵ0], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")], 
            dfs[!,Symbol("ΔSDS_Mx_x_0π")] .+ dfs[!,Symbol("ΔSDS_Mx_x_π0")],
            df_ϕs[!,Symbol("ΔSDS_A1_00")]*Ls^2; linewidth=er_lw, 
            whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:ϵ0], dfs[!,Symbol("SDS_Mx_x_0π")] .+ dfs[!,Symbol("SDS_Mx_x_π0")] ;
        marker = :xcross, color = colorschemes[:tab20][1-idx/num_points], 
        label = "β=$((1/Ts))", markersize=10)

    CairoMakie.scatter!(top, df_ϕs[!,:ϵ0], df_ϕs[!,Symbol("SDS_A1_00")]*Ls^2 , 
        marker = '□', color = colorschemes[:tab20][1-idx/num_points], 
         markersize=15)      
end

################
## Nematic order parameter as a function of strain ϵ₀
################
top=Axis(fig[3, 1], xticklabelsvisible=false, ylabel=L"$\left| \Phi^{B_1} \right|$", 
    ylabelsize=30, limits= ((-0.005, 0.25), nothing),  xlabelsize=30)

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    num_points=length(df_OP_LT[:,:T]);
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP_U0)) 
 
    eb=errorbars!(df_OPs[!,:ϵ0], df_OPs[!,Symbol("B1_OP")], 
        df_OPs[!,Symbol("ΔB1_OP")];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:ϵ0], df_OPs[!,Symbol("B1_OP")] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_points], label = "β=$((1/Ts))", markersize=10) 
         
end
########################
## Magnetic Order parameter, as a function of temperature
########################
top=Axis(fig[4, 1], xticklabelsvisible=false, ylabel=L"$\left| \left(  \mathbf{M}_{\mathbf{Q}_1}, \!
    \mathbf{M}_{\mathbf{Q}_2}  \right)  \right|$", 
    ylabelsize=30, limits= ((-0.005, 0.25), nothing),  xlabelsize=30)
CairoMakie.translate!(vlines!(top, [1.0/3, 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_OP_LT[:,:T],df_OP_LT[:,:L]))
    num_points=length(df_OP_LT[:,:T]);
    df_OPs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_OP_U0)) 

    eb=errorbars!(df_OPs[!,:ϵ0], df_OPs[!,Symbol("Mx_X_OP34")], df_OPs[!,Symbol("Δ" * "Mx_X_OP34")]; linewidth=0.5, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, df_OPs[!,:ϵ0], df_OPs[!,Symbol("Mx_X_OP34")] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_points], 
        label = "β=$((1/Ts))", markersize=10)  
end
################
## Superconducting pairing susceptibility as a function of strain ϵ₀
################
top=Axis(fig[5, 1], xlabel=L"$ϵ_0$", xlabelsize=30,
    ylabel=L"$χ_{\mathrm{pair}}^{zz} [\mathbf{Q}_0]$", 
    ylabelsize=30, limits= ((-0.005, 0.25), nothing))
    

CairoMakie.translate!(vlines!(top, [ 0.7, 1], color = :gray, linewidth=0.2), 0, 0, -0.8)
CairoMakie.translate!(hlines!(top, [0.2, 0.4], color = :gray,linewidth=0.2), 0, 0, -0.8)

for (idx, (Ts, Ls)) in enumerate(zip(df_LT[:,:T],df_LT[:,:L]))
    num_points=length(df_LT[:,:T]);
    dfs=filter(:L=>L -> L==Ls ,filter(:T =>T ->T==Ts, df_U0))  

    eb=errorbars!(dfs[!,:ϵ0], dfs[!,:PDS_spm_00], dfs[!,:ΔPDS_spm_00]; linewidth=er_lw, 
        whiskerwidth=10,color = colorschemes[:tab20][1-idx/num_points])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(top, dfs[!,:ϵ0], dfs[!,:PDS_spm_00] ;  marker = :xcross,
        color = colorschemes[:tab20][1-idx/num_points], 
        label = "β=$(1/Ts)", markersize=10)      
end




display(fig)
CairoMakie.save(joinpath(fig_path, "chi_vs_strain_U_" * to_string(U0) * ".png"), fig)

end