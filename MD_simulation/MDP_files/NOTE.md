# GROMACS parameter files

This directory contains reference input files used in the MD workflow.

- `em_real.mdp`: energy minimization
- `NVT.mdp`: 1 ns NVT equilibration
- `NPT.mdp`: 2 ns NPT equilibration
- `MD.mdp`: 200 ns production MD
- `ions.mdp`: ion-preparation input when required
- `mmpbsa.in`: MM-PBSA input
- `packmol.inp` and `topol.top`: a valid P2-C18 + mPEG2k-PCL2k pairwise reference example using files deposited under `../pairwise_prodrug_polymer_MD/training/P2-C18/`

The individual simulation directories contain the local molecular coordinates, topology/parameter files, `packmol.inp`, `topol.top` and `.mdp` files used for each system. The root batch scripts read those local files rather than the reference example in this directory.
