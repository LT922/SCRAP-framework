# Prodrug self-association analysis

Corresponding simulations:

```text
../../MD_simulation/prodrug_self_association_MD/
├── mechanistic_controls/
├── training/
└── validation/
```

Run one system at a time and explicitly provide the `mol` key from `segment_defination.tsv`.

P2-C18 example:

```bash
bash prepare_mmpbsa_index.sh \
  ../../MD_simulation/prodrug_self_association_MD/training/P2-C18 \
  ptx2c18

bash mmpbsa_analysis.sh \
  ../../MD_simulation/prodrug_self_association_MD/training/P2-C18 \
  ptx2c18
```

## Original simulation index

For a prodrug, the MD batch script creates:

```text
2       residue1
3       residue2
7       drug = residue1 + residue2
```

For a parent drug:

```text
2       parent-drug residues
6       drug
```

## Molecule-specific MM-PBSA index

The self-association box contains two copies of the same molecule.

`prepare_mmpbsa_index.sh` does not read any ITP file. Instead it uses GROMACS `splitres` on the existing residue groups.

For a prodrug, starting from groups 2 and 3:

```text
splitres 2  -> residue1_mol1, residue1_mol2
splitres 3  -> residue2_mol1, residue2_mol2
```

The script then combines the corresponding residue groups to create:

```text
Pro1
Pro2
```

`segment_defination.tsv` determines which split residue groups are Par1/Par2 or Mod1/Mod2.

The final four SCRAP self-association terms are:

```text
Pro-Pro
Par-Par
Mod-Mod
Par-Mod
```
