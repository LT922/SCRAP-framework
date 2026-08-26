# Structure preparation

This directory documents how the molecular inputs used in the MD simulations were generated.

## Small molecules

Drug and prodrug structures were prepared using:

```text
ChemDraw 22.0.0
→ Chem3D
→ MM2 minimization (minimum RMS gradient 0.001)
→ MOL2
→ CGenFF
→ cgenff_charmm2gmx.py 2.1.0
→ GROMACS-compatible structure and topology files
```

The `.mol2`, `.str`, `.pdb`, `.itp` and `.prm` files retained in the individual simulation directories document the resulting structures and force-field parameters.

## Polymer

The mPEG2k-PCL2k carrier required a separate reconstruction and equilibration procedure. The source template, ztop reconstruction procedure and final polymer input are documented in:

```text
polymer/NOTE.md
```
