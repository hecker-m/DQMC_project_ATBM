# IsingX_Neel
DQMC implementation based on `MonteCarlo.jl` package.
In a first little exercise, I reproduce results from https://doi.org/10.1103/PhysRevLett.125.247001.

Hereby, I slightly modify the `MonteCarlo.jl` package such that it also runs with continuous configuration fields.
For more details, please see the docs.
The physical problem considered in this project starts from the interaction Hamiltonian 

$ \hat{H}_{int} = -\frac{U}{N_ϕ}  \sum_{\boldsymbol{i}} \sum_{ζ=1}^{N_ϕ}  \, \hat{M}_{\boldsymbol{i}}^ζ \, \hat{M}_{\boldsymbol{i}}^ζ $,

similar as in Ref.[1,2]
