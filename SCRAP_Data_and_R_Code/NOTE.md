# Repository Notes

## Overview

This repository contains the source data, R scripts, molecular-dynamics (MD) analyses, processed tables, and figures used for the SCRAP prodrug–polymer compatibility study.

The repository is organized by analysis type rather than manuscript figure number.

```text
data/      Source data
R/         Analysis scripts
results/   Generated figures and tables
archive/   Historical or exploratory files
```

Main analysis domains:

```text
coassembly
combination
compusyn
md
training
validation
```

---

## 1. `data/`

Contains source data used by the active R workflows.

### `data/coassembly/`

- `combination.xlsx` — binary and multi-prodrug co-assembly data.

### `data/combination/`

- `combination_index.xlsx` — in vitro combination screening, IC50, viability, and CI data.
- `combination_therapy.xlsx` — tumour growth, body weight, and tumour-growth inhibition data.

### `data/compusyn/`

- `report.html` — original CompuSyn report containing dose–effect, median-effect, Fa–CI, CI, and DRI outputs.

### `data/md/`

#### `training_library_MD/`

- `MD_drug1polymer1.xlsx` — pairwise prodrug–polymer MD summary data.
- `MD_drugspolymers.xlsx` — multi-molecule MD summary data.

#### `DTX_Abi_coassembly_MD/`

Raw GROMACS analysis data for DTX-SI-C18, Abi-SI-C18, and the mixed DTX-SI-C18/Abi-SI-C18 system.

Common raw outputs:

- `rmsd*` — RMSD;
- `*gyrate*` — radius of gyration;
- `*-dist`, `*-cogdist`, `*-mindist` — distance metrics;
- `*-contacts` — contact number;
- `hbnum-*` — hydrogen bonds;
- `rdf-*` — radial distribution functions;
- `sasa-*`, `area_*` — SASA/contact-area analysis;
- `temperature`, `pressure`, `density`, `potential` — MD quality control;
- `size`, `ntr`, `nclid`, `cluster*` — clustering analysis.

Snapshot folders contain representative MD structures at selected time points.

### `data/publication/`

- `Supplementary_Data.xlsx` — publication-facing supplementary data workbook.

### `data/training/`

- `training_size_pdi.csv` — DLS size and PDI of the training library.
- `Training_segment_data.csv` — 13 segment-resolved SCRAP interaction descriptors.
- `bulk_descriptor_training_lib.csv` — whole-molecule descriptor data.

### `data/validation/`

- `test_lib.xlsx` — independent test-library data used for validation.

---

## 2. `R/`

Contains the active analysis scripts.

### `R/coassembly/`

- `Coassembly_formulation.R` — formulation-level co-assembly, rescue, size, and score analyses.

### `R/combination/`

- `Combination_therapy.R` — in vitro combination screening and in vivo combination-therapy analysis.

### `R/compusyn/`

- `CompuSyn_analysis.R` — parses `report.html` and reconstructs dose–effect, median-effect, Fa–CI, CI50, ED-level CI, and isobologram-style analyses.

### `R/md/`

- `Training_library_MD.R` — pairwise and multi-molecule MD analysis for representative training-library systems.
- `DTX_Abi_coassembly_MD.R` — main mechanistic MD analysis of DTX-SI-C18/Abi-SI-C18 co-assembly.
- `DTX_Abi_coassembly_MD_QC.R` — extended MD quality-control and trajectory-level analysis.

### `R/training/`

- `Training_library.R` — training-library size/PDI summary.
- `SCRAP_framework.R` — core 13-term SCRAP scoring, favourable windows, importance, and failure-mode analysis.
- `SCRAP_model_validation.R` — robustness analysis including ablation, LOOCV, threshold perturbation, correlation, and PCA.

### `R/validation/`

- `Independent_validation.R` — independent test-library validation and whole-molecule descriptor benchmarking.

---

## 3. `results/figures/`

Contains generated figures.

### `coassembly/`

Co-assembly and rescue plots, including DTX-fraction size heatmaps, rescue analyses, multi-prodrug co-assembly size, and formulation summaries.

### `combination/`

- `screening/` — IC50 curves, viability analyses, fixed-ratio screening, and CI comparisons.
- `therapy/` — tumour growth, tumour-growth inhibition, body weight, and treatment summaries.

### `compusyn/`

CompuSyn-derived dose–effect, median-effect, Fa–CI, CI50, ED-level CI, isobologram, and polygonogram plots.

### `md/training_library_MD/`

- `pairwise/` — RMSD, Rg, distance, hydrogen bonds, contact area, and carrier-integration summaries.
- `multi_molecule/` — RMSD, structural/energy summaries, hydrogen bonds, RDF, and metric heatmaps.

### `md/DTX_Abi_coassembly_MD/`

- `summary/` — compact mechanistic summary plots.
- `mechanistic_support/` — focused mechanistic analyses.
- `trajectory_support/` — trajectory-level support plots.
- `analysis/` — complete MD analysis outputs.
- `QC_support/` — independent QC and individual trajectory diagnostics.

### `training/`

- `training_library/` — DLS size/PDI overview.
- `SCRAP_framework/` — interaction-energy matrix, favourable windows, score, importance, and failure-mode plots.
- `SCRAP_model_validation/` — robustness, ablation, threshold perturbation, correlation, and PCA plots.

### `validation/independent_validation/`

Independent-validation plots including score versus DLS size, confusion matrix, window-satisfaction map, failed-window decomposition, module retention, and descriptor benchmarking.

---

## 4. `results/tables/`

Contains processed numerical outputs corresponding to the figures.

### `coassembly/`

Co-assembly size, rescue, and score-summary tables.

### `combination/`

- `screening/` — IC50, viability, fixed-ratio, and CI tables.
- `therapy/` — tumour growth, inhibition, body-weight, and treatment-summary tables.

### `compusyn/`

Parsed CompuSyn dose–effect, median-effect, Fa–CI, CI50, and ED-level CI data.

### `md/`

- `training_library_MD/` — pairwise and multi-molecule MD summary tables.
- `DTX_Abi_coassembly_MD/analysis/` — converted XVG files and processed MD summaries.
- `DTX_Abi_coassembly_MD/QC_support/` — QC-derived MD summaries.

### `training/`

- `training_library/` — processed size/PDI data.
- `SCRAP_framework/` — scores, failed windows, and importance tables.
- `SCRAP_model_validation/` — LOOCV, threshold, PCA, correlation, ablation, and performance tables.

### `validation/independent_validation/`

Test-library SCRAP scores, confusion matrices, module summaries, and whole-molecule descriptor benchmark tables.

---

## 5. `archive/`

Contains non-canonical historical or exploratory material.

- `exploratory_mixed_score/` — exploratory mixed-score analysis.
- `legacy_scripts/` — older plotting scripts.
- `pre_numberless_outputs_*` — backup of the previous figure-number-based code structure.

These files are not required for the active workflow.

---

## 6. SCRAP terminology

Canonical labels:

- `Pro` — intact prodrug;
- `Par` — parent-drug moiety;
- `Mod` — modifier/linker segment;
- `Pol` — whole polymer;
- `PEG` — PEG block;
- `PCL` — PCL block.

The 13 SCRAP interaction terms are:

```text
Pro–Pol
Par–Pol
Mod–Pol
Pro–PCL
Par–PCL
Mod–PCL
Pro–PEG
Par–PEG
Mod–PEG
Pro–Pro
Par–Par
Mod–Mod
Par–Mod
```

Raw MD filenames retain the original GROMACS naming where necessary.

---

## 7. Repository mapping

```text
data/coassembly/   → R/coassembly/   → results/{figures,tables}/coassembly/
data/combination/  → R/combination/  → results/{figures,tables}/combination/
data/compusyn/     → R/compusyn/     → results/{figures,tables}/compusyn/
data/md/           → R/md/           → results/{figures,tables}/md/
data/training/     → R/training/     → results/{figures,tables}/training/
data/validation/   → R/validation/   → results/{figures,tables}/validation/
```

Generated outputs use descriptive names instead of manuscript Figure or Supplementary Figure numbers.

