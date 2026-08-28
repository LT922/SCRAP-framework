# Multi-molecule analysis

Corresponding simulations:

```text
../../MD_simulation/multi_molecule_MD/
```

These simulations are used for assembly-level structural analysis.

## Standard single-species systems

Examples:

```bash
bash gromacs_analysis.sh \
  ../../MD_simulation/multi_molecule_MD/P2-C18 \
  ptx2c18

bash gromacs_analysis.sh \
  ../../MD_simulation/multi_molecule_MD/DTX-SI-C18 \
  dtxsic18

bash gromacs_analysis.sh \
  ../../MD_simulation/multi_molecule_MD/Abi-SI-C18 \
  abisic18
```

The standard index layout is identical to the pairwise MD layout:

Prodrug:

```text
7/8     residue1/residue2
12      Pol
13      Pro
14      Pol+Pro
15      PEG
16      PCL
```

Parent drug:

```text
7       Pro = Par
11      Pol
12      Pol+Pro
13      PEG
14      PCL
```

## Mixed DTX-SI-C18 + Abi-SI-C18

Use the dedicated script:

```bash
bash gromacs_analysis_DTX_Abi.sh \
  ../../MD_simulation/multi_molecule_MD/DTX-SI-C18+Abi-SI-C18 \
  dtxsic18 \
  abisic18
```

The current mixed-system MD index command gives:

```text
7/8      DTX-SI-C18 residue1/residue2
9/10     Abi-SI-C18 residue1/residue2
14       Pol
15       DTX-SI-C18
16       Abi-SI-C18
17       total Pro
18       Pol+Pro
19       PEG
20       PCL
```

According to `segment_definition.tsv`:

```text
dtxsic18  DTL  DTX  mod-par
abisic18  ABI  ABL  par-mod
```

Therefore:

```text
DTX-SI-C18: Mod=7, Par=8, Pro=15
Abi-SI-C18: Par=9, Mod=10, Pro=16
```

## Output compatibility

The multi-molecule scripts retain the historical XVG stems used by the deposited R workflow. In particular, the mixed DTX-SI-C18/Abi-SI-C18 script regenerates:

```text
whole-component COM and COG distances
segment-resolved COM distances
whole- and segment-resolved H-bonds
DTX/Abi/polymer/PCL/PEG contact outputs
all seven deposited RDF outputs
all deposited SASA component/union outputs
drug-, DTX- and Abi-specific cluster outputs
```

This keeps `Post_MD_analysis` consistent with the raw-file dictionaries used by the DTX/Abi coassembly R analysis.

## RDF settings

RDFs are calculated with molecular centres of mass for both reference and selected groups, matching the deposited SI analysis:

```text
-selrpos mol_com
-seltype mol_com
-pbc yes
-bin 0.002
```

Both RDF (`-o`) and cumulative RDF (`-cn`) outputs are generated.
