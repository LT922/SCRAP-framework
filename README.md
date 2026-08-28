# Reproducibility repository

This repository contains the computational workflow and analysis code for the segment-resolved compatibility framework (SCRAP). It is organized into three sequential parts: **MD simulation**, **post-MD analysis**, and **data analysis/figure generation**.

## Repository structure

```text
.
├── MD_simulation/
├── Post_MD_analysis/
└── SCRAP_Data_and_R_Code/
```

### 1. `MD_simulation/`

Contains the molecular-dynamics input framework, including:

- molecular and polymer structure preparation;
- CHARMM36/CGenFF-compatible force-field files;
- GROMACS parameter files;
- pairwise prodrug–polymer simulations;
- prodrug self-association simulations;
- multi-molecule assembly simulations;
- batch scripts for the three simulation classes.

The corresponding `NOTE.md` files describe structure preparation, polymer construction, simulation settings and execution.

### 2. `Post_MD_analysis/`

Contains the post-processing workflow corresponding directly to the three MD simulation classes.

```text
MD_simulation/
        ↓
Post_MD_analysis/
```

GROMACS scripts are provided for structural analyses, and `gmx_MMPBSA` scripts are provided for the segment-resolved interaction-energy calculations used in SCRAP.

`segment_definition.tsv` contains the manually verified assignment of prodrug residues to the parent-drug (**Par**) and modifier/linker (**Mod**) segments. Large regenerated trajectories and intermediate analysis files are written to `analysis_output/` and are not intended for version control.

### 3. `SCRAP_Data_and_R_Code/`

Contains the compact datasets and R scripts used for SCRAP scoring, validation, statistical analysis and figure generation.

```text
Post_MD_analysis/
        ↓
compact analysis data
        ↓
R scripts
        ↓
figures and tables
```

The `data/` directory contains deposited analysis inputs, `R/` contains the corresponding analysis scripts, and regenerated outputs are written under `results/`.

## Workflow

```text
Structure preparation
        ↓
MD simulation
        ↓
GROMACS post-processing
        ↓
segment-resolved MM-PBSA
        ↓
SCRAP scoring and validation
        ↓
statistical analysis and figure generation
```

The segment notation used throughout the repository is:

```text
Pro  intact prodrug
Par  parent-drug moiety
Mod  modifier/linker moiety
Pol  complete mPEG2k-PCL2k polymer
PEG  PEG block
PCL  PCL block
```

Large raw MD trajectories and regenerated intermediate outputs are intentionally excluded from the repository. The deposited files provide the simulation inputs, post-processing commands, compact analysis data and R code required to reproduce the computational workflow. Detailed instructions are provided in the `NOTE.md` files within each subdirectory.
