# DQMC_project_ATBM
A DQMC project based on the (attractive two-band model---ATBM) Hamiltonian introduced in https://doi.org/10.1103/PhysRevLett.125.247001. First, in `IsingX_Neel` I reproduce parts of the paper's results. Later, in `IsingX_Striped` and `Heisenberg_Striped` I choose a Fermi surface which favors a degenerate magnetic ground state with $N_\phi=1$ and $N_\phi=3$ bosonic components. These cases are studied with respect to all sorts of bilinear orders.

The implementation is based on the `MonteCarlo.jl` package, which I modified slightly. 

Research in progress.
