# Pairwise prodrug-polymer analysis

Corresponding simulations:

```text
../../MD_simulation/pairwise_prodrug_polymer_MD/
├── mechanistic_controls/
├── training/
└── validation/
```

Run one system at a time and explicitly provide its `mol` key from `../segment_defination.tsv`.

P2-C18 example:

```bash
bash gromacs_analysis.sh \
  ../../MD_simulation/pairwise_prodrug_polymer_MD/training/P2-C18 \
  ptx2c18

bash mmpbsa_analysis.sh \
  ../../MD_simulation/pairwise_prodrug_polymer_MD/training/P2-C18 \
  ptx2c18
```

## Index layout

For a prodrug, the MD batch script creates:

```text
7       residue1
8       residue2
12      polymer
13      drug
14      polymer&drug
15      peg
16      pcl
```

The TSV decides whether group 7/8 are Par/Mod.

For a parent drug, the MD batch script creates:

```text
7       parent drug
11      polymer
12      polymer&drug
13      peg
14      pcl
```

Therefore:

```text
Pro = Par = 7
Pol = 11
PEG = 13
PCL = 14
```

## GROMACS post-processing

`gromacs_analysis.sh` performs:

```text
PBC correction
RMSD
Rg
COM/COG distances
hydrogen bonds
SASA
clustering
representative trajectory
```

Commands are written explicitly in the script so the selected groups and output files can be read directly.

## MM-PBSA

`mmpbsa_analysis.sh` calculates the nine carrier-facing SCRAP terms for a prodrug:

```text
Pro-Pol   Par-Pol   Mod-Pol
Pro-PCL   Par-PCL   Mod-PCL
Pro-PEG   Par-PEG   Mod-PEG
```

For a parent-drug control, only whole-parent interactions are calculated.

## MM-PBSA group order

To reproduce the original calculation convention, the `-cg` order is:

```text
carrier first, drug/prodrug segment second
```

For example, P2-C18 Par-PCL is run as:

```bash
-cg 16 8
```

not `-cg 8 16`. This follows the original MM-PBSA batch workflow.
