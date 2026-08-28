# Pairwise prodrug-polymer MD

This directory contains the pairwise systems used to evaluate interactions between one drug/prodrug molecule and one mPEG2k-PCL2k chain.

```text
System composition: 1 drug/prodrug + 1 mPEG2k-PCL2k
Box:                8 × 8 × 8 nm³
Packmol output:     polymerdrug.pdb
```

Production duration: **200 ns** for every deposited production MD system in this class.

The polymer coordinate supplied to Packmol is the pre-equilibrated:

```text
mpeg2kpcl2knpt4ns.pdb
```

which was obtained from the ztop-built mPEG2k-PCL2k structure after GROMACS energy minimization, 2 ns NVT and 4 ns NPT equilibration.

Subdirectories:

- `training/`: systems used to establish the compatibility framework.
- `validation/`: independent validation systems.
- `mechanistic_controls/`: parent-drug/control systems used for mechanistic comparison.

Run from the repository root with:

```bash
bash MD_simulation/run_pairwise_prodrug_polymer_MD.sh
```

The batch script builds the system with Packmol, creates an 8 nm cubic box, solvates the system, generates the required index groups, and runs EM → NVT → NPT → production MD.
