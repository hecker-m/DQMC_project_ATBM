# DQMC_project_ATBM
A DQMC project based on the (attractive two-band model---ATBM) Hamiltonian introduced in https://doi.org/10.1103/PhysRevLett.125.247001, akin to an attractive two-band Hubbard model. First, in `IsingX_Neel` I reproduce parts of the paper's results. Second, in `IsingX_Striped` and `Heisenberg_Striped` (corresponding to $N_\phi=1$ and $N_\phi=3$ bosonic components, respectively) I choose a Fermi surface which favors the degenerate ordering vectors $\boldsymbol{Q}_1=(\pi, 0)$ and $\boldsymbol{Q}_2=(0, \pi)$.
Due to this degeneracy, bilinear orders such as nematicity are possible. 
These systems are studied with respect to both, primary order (magnetism) and all associated bilinear orders (nematic, translational-symmetry breaking, ...)

The implementation is based on the `MonteCarlo.jl` package. 

Peer-reviewed publication soon to be expected.
