# GROMACS parameter files

This directory contains reference GROMACS input files used in the MD workflow.

- `em_real.mdp`: energy minimization
- `NVT.mdp`: NVT equilibration
- `NPT.mdp`: NPT equilibration
- `MD.mdp`: production MD
- `ions.mdp`: ion-preparation input when required
- `mmpbsa.in`: MM-PBSA input
- `packmol.inp` and `topol.top`: reference templates

The individual simulation directories contain the local copies used for each system. The root batch scripts read the local files in each simulation directory rather than these reference templates directly.
