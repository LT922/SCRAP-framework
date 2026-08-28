# Post-MD analysis

This directory contains the post-processing workflow corresponding directly to the three MD simulation classes in `../MD_simulation/`.

```text
Post_MD_analysis/
├── NOTE.md
├── segment_definition.tsv
├── pairwise_prodrug_polymer_MD/
│   ├── NOTE.md
│   ├── gromacs_analysis.sh
│   └── mmpbsa_analysis.sh
├── prodrug_self_association_MD/
│   ├── NOTE.md
│   ├── prepare_mmpbsa_index.sh
│   └── mmpbsa_analysis.sh
├── multi_molecule_MD/
│   ├── NOTE.md
│   ├── gromacs_analysis.sh
│   └── gromacs_analysis_DTX_Abi.sh
└── analysis_output/
    └── NOTE.md
```

The scripts do **not** read the small-molecule ITP files. The manually verified `segment_definition.tsv` is the only source used to decide whether residue1/residue2 are Par or Mod.

## Segment definitions

The SCRAP terminology is:

```text
Pro = intact prodrug
Par = parent-drug moiety
Mod = modifier/linker moiety
Pol = intact mPEG2k-PCL2k
PEG = PEG block
PCL = PCL block
```

For each analysis, the user supplies:

```text
1. the MD system directory
2. the `mol` key from segment_definition.tsv
```

Example:

```bash
bash pairwise_prodrug_polymer_MD/gromacs_analysis.sh \
  ../MD_simulation/pairwise_prodrug_polymer_MD/training/P2-C18 \
  ptx2c18
```

The TSV row is:

```text
ptx2c18    PTL    PTX    mod-par
```

Therefore:

```text
residue1 = PTL = Mod
residue2 = PTX = Par
```

The current pairwise/prodrug MD index command defines:

```text
7  = residue1
8  = residue2
12 = Pol
13 = Pro
14 = Pol+Pro
15 = PEG
16 = PCL
```

Thus P2-C18 is analysed as:

```text
Mod = 7
Par = 8
Pro = 13
Pol = 12
PEG = 15
PCL = 16
```

The scripts print this mapping before any analysis starts and also save it to `group_mapping.log`.

## Relationship to SCRAP

Pairwise prodrug-polymer MD provides nine carrier-facing MM-PBSA terms:

```text
Pro-Pol   Par-Pol   Mod-Pol
Pro-PCL   Par-PCL   Mod-PCL
Pro-PEG   Par-PEG   Mod-PEG
```

Self-association MD provides four cohesion terms:

```text
Pro-Pro   Par-Par   Mod-Mod   Par-Mod
```

For prodrugs, `Par-Mod` is defined as the arithmetic mean of the reciprocal cross-molecule interactions `Par1-Mod2` and `Par2-Mod1`.

Together these are the 13 SCRAP descriptors.

Multi-molecule MD is used for assembly-level structural analysis and mechanistic support rather than for generating the 13-term score.

## Output

Generated trajectories and analysis files are written to:

```text
Post_MD_analysis/analysis_output/
```

These outputs can be regenerated from local MD trajectories and are not intended to be stored in GitHub.

## Software

The scripts correspond to:

- GROMACS 2023
- gmx_MMPBSA 1.6.4
- MPI-enabled gmx_MMPBSA

MM-PBSA uses:

```text
../MD_simulation/MDP_files/mmpbsa.in
```


