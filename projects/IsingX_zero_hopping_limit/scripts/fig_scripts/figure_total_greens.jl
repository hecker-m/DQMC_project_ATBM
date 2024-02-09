
#   fieldnames(typeof(mc.model.l))
N=length(lattice(mc))
μ0=2.0;
U0=mc.model.U
Ls=mc.model.l.Ls[1]
δτ=mc.parameters.delta_tau
β0=mc.parameters.beta
full_greens=mean(mc.measurements[:total_greens].observable)
full_greens_error=std_error(mc.measurements[:total_greens].observable)
M=size(full_greens,1)

G0ls=zeros(M, N, 4, 2)
for τ in 1:M
    for i in 1:N
        G0ls[τ, i, 1, 1]=full_greens[τ, 2, i, i]
        G0ls[τ, i, 2, 1]=full_greens[τ, 2, i, i+N]
        G0ls[τ, i, 3, 1]=full_greens[τ, 2, i+N, i]
        G0ls[τ, i, 4, 1]=full_greens[τ, 2, i+N, i+N]
        G0ls[τ, i, 1, 2]=full_greens_error[τ, 2, i, i]
        G0ls[τ, i, 2, 2]=full_greens_error[τ, 2, i, i+N]
        G0ls[τ, i, 3, 2]=full_greens_error[τ, 2, i+N, i]
        G0ls[τ, i, 4, 2]=full_greens_error[τ, 2, i+N, i+N]
    end
end
Gl0s=zeros(M, N, 4, 2)
for τ in 1:M
    for i in 1:N
        Gl0s[τ, i, 1, 1]=full_greens[τ, 3, i, i]
        Gl0s[τ, i, 2, 1]=full_greens[τ, 3, i, i+N]
        Gl0s[τ, i, 3, 1]=full_greens[τ, 3, i+N, i]
        Gl0s[τ, i, 4, 1]=full_greens[τ, 3, i+N, i+N]
        Gl0s[τ, i, 1, 2]=full_greens_error[τ, 3, i, i]
        Gl0s[τ, i, 2, 2]=full_greens_error[τ, 3, i, i+N]
        Gl0s[τ, i, 3, 2]=full_greens_error[τ, 3, i+N, i]
        Gl0s[τ, i, 4, 2]=full_greens_error[τ, 3, i+N, i+N]
    end
end
Glls=zeros(M, N, 4, 2)
for τ in 1:M
    for i in 1:N
        Glls[τ, i, 1, 1]=full_greens[τ, 4, i, i]
        Glls[τ, i, 2, 1]=full_greens[τ, 4, i, i+N]
        Glls[τ, i, 3, 1]=full_greens[τ, 4, i+N, i]
        Glls[τ, i, 4, 1]=full_greens[τ, 4, i+N, i+N]
        Glls[τ, i, 1, 2]=full_greens_error[τ, 4, i, i]
        Glls[τ, i, 2, 2]=full_greens_error[τ, 4, i, i+N]
        Glls[τ, i, 3, 2]=full_greens_error[τ, 4, i+N, i]
        Glls[τ, i, 4, 2]=full_greens_error[τ, 4, i+N, i+N]
    end
end

path_Analysis_Fcns="/home/mhecker/Google Drive/DQMC/AFM_2_band_model/DQMC_project_ATBM/src/Analysis_Fcns/"
include(path_Analysis_Fcns * "analytical_fcns.jl")

using Plots, DataFrames, StatsPlots, Query, FileIO, CairoMakie, ColorSchemes, LaTeXStrings, CSV
er_lw=0.2;
#=
################
## G_{ii}(0,τ)
################
for i0 in [1,2,16]
fig = Figure(resolution = (800, 800))

f1=Axis(fig[1, 1], ylabel=L"$G_{11}(0,τ)$", title=L"\mathrm{position: }\;\; (%$(i0),%$(i0))",
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
f2=Axis(fig[2, 1], ylabel=L"$G_{12}(0,τ)$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
f3=Axis(fig[3, 1], ylabel=L"$G_{21}(0,τ)$", 
    xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
f4=Axis(fig[4, 1], xlabel=L"$τ*t$", ylabel=L"$G_{22}(0,τ)$", 
    xlabelsize=30, ylabelsize=30, 
    xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))  
rangeτ=range(0, β0, 30)

count=1;

l0=1;
eb=errorbars!(f1, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M],  [G0ls[τ, i0, l0, 2] for τ in 1:M];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
CairoMakie.translate!(eb, 0, 0, -0.5)
CairoMakie.scatterlines!(f1, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M];  
    marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
l0=2;
eb=errorbars!(f2, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M],  [G0ls[τ, i0, l0, 2] for τ in 1:M];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
CairoMakie.translate!(eb, 0, 0, -0.5)
CairoMakie.scatterlines!(f2, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M];  
    marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
l0=3;
eb=errorbars!(f3, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M],  [G0ls[τ, i0, l0, 2] for τ in 1:M];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
CairoMakie.translate!(eb, 0, 0, -0.5)
CairoMakie.scatterlines!(f3, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M];  
    marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
l0=4;
eb=errorbars!(f4, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M],  [G0ls[τ, i0, l0, 2] for τ in 1:M];
        linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
CairoMakie.translate!(eb, 0, 0, -0.5)
CairoMakie.scatterlines!(f4, [(τ-1)*δτ for τ in 1:M], [G0ls[τ, i0, l0, 1] for τ in 1:M];  
    marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")


count+=2;
CairoMakie.scatterlines!(f1, rangeτ, G0τ_red_1_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
    marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
    color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
CairoMakie.scatterlines!(f2, rangeτ, G0τ_red_2_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
    marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
    color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
CairoMakie.scatterlines!(f3, rangeτ, G0τ_red_3_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
    marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
    color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
CairoMakie.scatterlines!(f4, rangeτ, G0τ_red_4_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
    marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
    color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
# axislegend(top, position=(0, 1))
# axislegend(mid, position=(0.2, 0))

display(fig)
end

################
## G_{ii}(τ,0)
################
for i0 in [1,2,16]
    fig = Figure(resolution = (800, 800))
    
    f1=Axis(fig[1, 1], ylabel=L"$G_{11}(τ,0)$", title=L"\mathrm{position: }\;\; (%$(i0),%$(i0))",
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f2=Axis(fig[2, 1], ylabel=L"$G_{12}(τ,0)$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f3=Axis(fig[3, 1], ylabel=L"$G_{21}(τ,0)$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f4=Axis(fig[4, 1], xlabel=L"$τ*t$", ylabel=L"$G_{22}(τ,0)$", 
        xlabelsize=30, ylabelsize=30, 
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))  
    rangeτ=range(0, β0, 30)
    
    count=1;
    
    l0=1;
    eb=errorbars!(f1, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M],  [Gl0s[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f1, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    l0=2;
    eb=errorbars!(f2, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M],  [Gl0s[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f2, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    l0=3;
    eb=errorbars!(f3, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M],  [Gl0s[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f3, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    l0=4;
    eb=errorbars!(f4, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M],  [Gl0s[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f4, [(τ-1)*δτ for τ in 1:M], [Gl0s[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    
    
    count+=2;
    CairoMakie.scatterlines!(f1, rangeτ, Gτ0_red_1_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    CairoMakie.scatterlines!(f2, rangeτ, Gτ0_red_2_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    CairoMakie.scatterlines!(f3, rangeτ, Gτ0_red_3_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    CairoMakie.scatterlines!(f4, rangeτ, Gτ0_red_4_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    
    display(fig)
    end

################
## G_{ii}(τ,τ)
################
for i0 in [1,2,16]
    fig = Figure(resolution = (800, 800))
    
    f1=Axis(fig[1, 1], ylabel=L"$G_{11}(τ,τ)$", title=L"\mathrm{position: }\;\; (%$(i0),%$(i0))",
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f2=Axis(fig[2, 1], ylabel=L"$G_{12}(τ,τ)$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f3=Axis(fig[3, 1], ylabel=L"$G_{21}(τ,τ)$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f4=Axis(fig[4, 1], xlabel=L"$τ*t$", ylabel=L"$G_{22}(τ,τ)$", 
        xlabelsize=30, ylabelsize=30, 
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))  
    rangeτ=range(0, β0, 30)
    
    count=1;
    
    l0=1;
    eb=errorbars!(f1, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M],  [Glls[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f1, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    l0=2;
    eb=errorbars!(f2, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M],  [Glls[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f2, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    l0=3;
    eb=errorbars!(f3, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M],  [Glls[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f3, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    l0=4;
    eb=errorbars!(f4, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M],  [Glls[τ, i0, l0, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f4, [(τ-1)*δτ for τ in 1:M], [Glls[τ, i0, l0, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")
    
    
    count+=2;
    CairoMakie.scatterlines!(f1, rangeτ, [Gτ0_red_1_analytic(U0, β0, 0.5*(μ0-(-μ0)), 0) for τ in rangeτ] ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    CairoMakie.scatterlines!(f2, rangeτ, [Gτ0_red_2_analytic(U0, β0, 0.5*(μ0-(-μ0)), 0) for τ in rangeτ] ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    CairoMakie.scatterlines!(f3, rangeτ, [Gτ0_red_3_analytic(U0, β0, 0.5*(μ0-(-μ0)), 0) for τ in rangeτ] ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    CairoMakie.scatterlines!(f4, rangeτ, [Gτ0_red_4_analytic(U0, β0, 0.5*(μ0-(-μ0)), 0) for τ in rangeτ] ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    
    display(fig)
    end

=#

################
## G_{ii}(τ,0)
################
for i0 in [1,], j0 in [1,3]
    fig = Figure(resolution = (800, 800))
    
    f1=Axis(fig[1, 1], ylabel=L"$G_{11}^{%$(i0),%$(i0)}(τ,0)$", title=L"\mathrm{position: }\;\; (%$(i0),%$(i0))",
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f2=Axis(fig[2, 1], ylabel=L"$G_{11}^{%$(j0),%$(j0)}(0,τ)$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
 
    rangeτ=range(0, β0, 30)
    
    count=1;
    valsτ=[(τ-1)*δτ for τ in 1:M]
    eb=errorbars!(f1, valsτ, [Gl0s[τ, i0, 1, 1] for τ in 1:M],  [Gl0s[τ, i0, 1, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f1, valsτ, [Gl0s[τ, i0, 1, 1] for τ in 1:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")

    eb=errorbars!(f2, valsτ[2:end], [G0ls[τ, j0, 1, 1] for τ in 2:M],  [G0ls[τ, j0, 1, 2] for τ in 2:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f2, valsτ[2:end], [G0ls[τ, j0, 1, 1] for τ in 2:M];  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")

    
    count+=2;
    CairoMakie.scatterlines!(f1, rangeτ, Gτ0_red_1_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    CairoMakie.scatterlines!(f2, rangeτ, G0τ_red_1_analytic.(U0, β0, 0.5*(μ0-(-μ0)), rangeτ) ;  
        marker =:xcross, markersize=1, linestyle=:dot, linewidth=2.5,
        color = colorschemes[:tab20][count], label = L"$(ana):\;β=%$(β0), L=%$(Ls) $")
    
    display(fig)
    end

################
## G_{ii}(τ,0)
################
for i0 in [1,2, 3], j0 in [1,2, 3]
    fig = Figure(resolution = (800, 800))
    
    f1=Axis(fig[1, 1], ylabel=L"$G_{11}^{%$(i0),%$(i0)}(τ,0)$", title=L"\mathrm{position: }\;\; (%$(i0),%$(i0))",
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
    f2=Axis(fig[2, 1], ylabel=L"$G_{11}^{%$(j0),%$(j0)}(0,τ)$", 
        xlabelsize=30, ylabelsize=30, xticklabelsvisible=false,
        xticklabelsize=20, yticklabelsize=20, limits= (nothing, nothing))
 
    rangeτ=range(0, β0, 30)
    CairoMakie.lines!(f1, rangeτ, zeros(length(rangeτ)) ; linewidth=1.5, color = :black) 
    CairoMakie.lines!(f2, rangeτ, zeros(length(rangeτ)) ; linewidth=1.5, color = :black) 

    count=1;
    valsτ=[(τ-1)*δτ for τ in 1:M]
    eb=errorbars!(f1, valsτ, [Gl0s[τ, i0, 1, 1] for τ in 1:M] .- Gτ0_red_1_analytic.(U0, β0, μ0, valsτ),  
            [Gl0s[τ, i0, 1, 2] for τ in 1:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f1, valsτ, [Gl0s[τ, i0, 1, 1] for τ in 1:M] .- Gτ0_red_1_analytic.(U0, β0, μ0, valsτ);  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")

    eb=errorbars!(f2, valsτ[2:end], [G0ls[τ, j0, 1, 1] for τ in 2:M] .- G0τ_red_1_analytic.(U0, β0, μ0, valsτ[2:end]),
            [G0ls[τ, j0, 1, 2] for τ in 2:M];
            linewidth=er_lw, whiskerwidth=10,color = colorschemes[:tab20][count])
    CairoMakie.translate!(eb, 0, 0, -0.5)
    CairoMakie.scatterlines!(f2, valsτ[2:end], [G0ls[τ, j0, 1, 1] for τ in 2:M] .- G0τ_red_1_analytic.(U0, β0, μ0, valsτ[2:end]);  
        marker =:xcross, markersize=7, color = colorschemes[:tab20][count], label = L"$():\;β=%$(β0), L=%$(Ls) $")

    
    display(fig)
    end
