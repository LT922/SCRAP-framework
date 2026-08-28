# Molecular dynamics input files and reproduction guide

This repository contains the molecular structures, force-field parameters and GROMACS input files used for the MD calculations.

## Repository layout

```text
.
├── MDP_files/
├── forcefields/
│   └── charmm36-jul2022.ff/
├── structure_preparation/
│   └── polymer/
│       ├── EG6CL6/
│       └── mpeg2kpcl2k/
├── pairwise_prodrug_polymer_MD/
├── prodrug_self_association_MD/
├── multi_molecule_MD/
├── run_pairwise_prodrug_polymer_MD.sh
├── run_prodrug_self_association_MD.sh
└── run_multi_molecule_MD.sh
```

- `structure_preparation/`: source structures and polymer-construction files.
- `MDP_files/`: reference GROMACS parameter files.
- `forcefields/`: shared CHARMM36 July 2022 force field.
- `pairwise_prodrug_polymer_MD/`: one drug/prodrug + one mPEG2k-PCL2k chain.
- `prodrug_self_association_MD/`: two copies of the same drug/prodrug.
- `multi_molecule_MD/`: multi-molecule assembly simulations.
- `run_*_MD.sh`: root-level batch scripts used to build and run all systems in each simulation class.

Each simulation directory contains the molecular coordinates, system-specific topology/parameter files, `packmol.inp`, `topol.top` and local `.mdp` files required to reproduce that system.

## Small-molecule preparation

The actual structure-preparation workflow used for drugs and prodrugs was:

```text
ChemDraw
→ SMILES
→ Chem3D
→ MM2 minimization
→ MOL2
→ CGenFF
→ cgenff_charmm2gmx.py
→ GROMACS-compatible PDB / ITP / PRM / TOP
```

Structures were drawn in **ChemDraw 22.0.0**. SMILES strings were converted to three-dimensional structures in Chem3D and minimized with the MM2 force field to a minimum RMS gradient of 0.001. The minimized structures were saved as `.mol2`, submitted to the CGenFF server, and the resulting `.str` files were converted to GROMACS format using `cgenff_charmm2gmx.py` (version 2.1.0).

No xTB geometry-optimization step was used for the structures in this repository.

### Historical local topology provenance

Each simulation directory preserves the local topology files used by that deposited system. Some self-association directories contain locally preprocessed topology files whose formatting or atom naming differs from the corresponding pairwise system. These local files are intentionally preserved rather than retrospectively replacing them with pairwise topologies, because such replacement would change the simulation input rather than merely correct repository organization.

## mPEG2k-PCL2k preparation

mPEG2k-PCL2k was reconstructed from the lower-molecular-weight `EG6CL6` template using the **ztop** Python package. Details and the exact ztop commands are given in:

```text
structure_preparation/polymer/NOTE.md
```

The ztop-generated full polymer was subsequently equilibrated in GROMACS:

```text
ztop-built mPEG2k-PCL2k
→ energy minimization
→ 2 ns NVT
→ 4 ns NPT
→ mpeg2kpcl2knpt4ns.pdb
```

**`mpeg2kpcl2knpt4ns.pdb` is the polymer coordinate file actually used to construct the pairwise and multi-molecule MD systems.** The corresponding polymer topology and parameters are supplied as `mpeg2kpcl2k.itp` and `mpeg2kpcl2k.prm`.

## MD software and hardware

Original calculations were performed using:

| Component | Version / specification |
|---|---|
| GROMACS | 2023 |
| Packmol | 20.11.0 |
| CHARMM force field | CHARMM36, July 2022 |
| cgenff_charmm2gmx.py | 2.1.0 |
| CUDA | 11.7 |
| CPU | Intel Core i9-12900K |
| GPU | NVIDIA GeForce RTX 3060 Ti |

GROMACS production runs used GPU acceleration for non-bonded calculations.

## Running the simulations

Run the batch scripts from the repository root (the scripts resolve their own location, so they do not depend on the current working directory):

```bash
bash MD_simulation/run_pairwise_prodrug_polymer_MD.sh
bash MD_simulation/run_prodrug_self_association_MD.sh
bash MD_simulation/run_multi_molecule_MD.sh
```

The scripts perform:

```text
Packmol
→ editconf
→ solvation
→ index generation
→ energy minimization
→ NVT
→ NPT
→ production MD
```

The scripts are restart-aware: completed systems are skipped and interrupted production runs can continue from `MD.cpt`.

See the local `NOTE.md` files inside each simulation-class directory for system composition and box size.
