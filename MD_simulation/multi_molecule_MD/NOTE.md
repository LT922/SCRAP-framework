# Multi-molecule MD

This directory contains assembly-level drug/prodrug-polymer simulations.

For the standard single-prodrug systems:

```text
System composition: 20 drug/prodrug molecules + 30 mPEG2k-PCL2k chains
Box:                14 × 14 × 14 nm³
Packmol output:     polymerdrug.pdb
```

For the mixed `DTX-SI-C18+Abi-SI-C18` system:

```text
10 DTX-SI-C18 + 10 Abi-SI-C18 + 30 mPEG2k-PCL2k chains
```

The polymer coordinate supplied to Packmol is:

```text
mpeg2kpcl2knpt4ns.pdb
```

This structure was obtained from the ztop-built mPEG2k-PCL2k chain after GROMACS energy minimization, 2 ns NVT and 4 ns NPT equilibration.

Run all multi-molecule systems from the repository root with:

```bash
./run_multi_molecule_MD.sh
```

The batch script uses a 14 nm cubic box and automatically applies the appropriate index-generation scheme for parent drugs, single-prodrug systems and the mixed DTX-SI-C18/Abi-SI-C18 system.
