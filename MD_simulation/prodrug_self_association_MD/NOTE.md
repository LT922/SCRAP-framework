# Prodrug self-association MD

This directory contains two-molecule systems used to evaluate drug/prodrug self-association.

```text
System composition: 2 copies of the same drug/prodrug
Box:                6 × 6 × 6 nm³
Packmol output:     drug2.pdb
```

Production duration: **200 ns** for every deposited production MD system in this class.

The directory is divided into:

- `training/`
- `validation/`
- `mechanistic_controls/`

No polymer is present in these simulations.

Run from the repository root with:

```bash
bash MD_simulation/run_prodrug_self_association_MD.sh
```

The batch script performs Packmol → box generation → solvation → index generation → EM → NVT → NPT → production MD.


## Local topology provenance

Use each system's deposited local `topol.top` and molecule topology files when reproducing that self-association simulation. Some local self-association topology files differ in formatting or atom naming from the corresponding pairwise-system files. They are retained as deposited simulation inputs and should not be replaced solely to make the two simulation classes textually identical.
