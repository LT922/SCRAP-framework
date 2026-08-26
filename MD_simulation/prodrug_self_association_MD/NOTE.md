# Prodrug self-association MD

This directory contains two-molecule systems used to evaluate drug/prodrug self-association.

```text
System composition: 2 copies of the same drug/prodrug
Box:                6 × 6 × 6 nm³
Packmol output:     drug2.pdb
```

The directory is divided into:

- `training/`
- `validation/`
- `mechanistic_controls/`

No polymer is present in these simulations.

Run all systems from the repository root with:

```bash
./run_prodrug_self_association_MD.sh
```

The batch script performs Packmol → box generation → solvation → index generation → EM → NVT → NPT → production MD.
