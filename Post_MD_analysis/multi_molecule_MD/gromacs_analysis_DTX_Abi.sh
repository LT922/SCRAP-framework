#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Mixed DTX-SI-C18 + Abi-SI-C18 multi-molecule analysis
#
# Usage:
#   bash gromacs_analysis_DTX_Abi.sh \
#     PATH_TO_DTX-SI-C18+Abi-SI-C18 \
#     dtxsic18 \
#     abisic18
#
# This script intentionally reproduces the historical XVG names
# used by the deposited downstream R analysis.
# ============================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 PATH_TO_SYSTEM DTX_MOL_KEY ABI_MOL_KEY"
    exit 1
fi

GMX="${GMX:-gmx}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$POST_ROOT/.." && pwd)"

SEGMENT_TABLE="$POST_ROOT/segment_definition.tsv"
OUTPUT_ROOT="$POST_ROOT/analysis_output/multi_molecule_MD"

SYSTEM_DIR="$(cd "$1" && pwd)"
DTX_KEY="$2"
ABI_KEY="$3"

if [ "$(basename "$SYSTEM_DIR")" != "DTX-SI-C18+Abi-SI-C18" ]; then
    echo "ERROR: this script is only for DTX-SI-C18+Abi-SI-C18." >&2
    exit 1
fi

OUT="$OUTPUT_ROOT/DTX-SI-C18+Abi-SI-C18"
mkdir -p "$OUT"

for FILE in MD.tpr MD.xtc index.ndx; do
    if [ ! -f "$SYSTEM_DIR/$FILE" ]; then
        echo "ERROR: missing $SYSTEM_DIR/$FILE" >&2
        exit 1
    fi
done

# ------------------------------------------------------------
# Read the two manually verified TSV rows.
# No ITP file is read.
# ------------------------------------------------------------
read_segment_row() {
    local key="$1"

    awk -F '\t' -v key="$key" '
        NR > 1 && $1 == key {
            print $2 "|" $3 "|" $4
        }
    ' "$SEGMENT_TABLE"
}

DTX_ROW="$(read_segment_row "$DTX_KEY")"
ABI_ROW="$(read_segment_row "$ABI_KEY")"

[ -n "$DTX_ROW" ] || { echo "ERROR: '$DTX_KEY' not found in TSV." >&2; exit 1; }
[ -n "$ABI_ROW" ] || { echo "ERROR: '$ABI_KEY' not found in TSV." >&2; exit 1; }

IFS='|' read -r DTX_RES1 DTX_RES2 DTX_ORDER <<< "$DTX_ROW"
IFS='|' read -r ABI_RES1 ABI_RES2 ABI_ORDER <<< "$ABI_ROW"

if [ "$DTX_ORDER" != "mod-par" ]; then
    echo "ERROR: expected DTX-SI-C18 to be mod-par, found '$DTX_ORDER'." >&2
    exit 1
fi

if [ "$ABI_ORDER" != "par-mod" ]; then
    echo "ERROR: expected Abi-SI-C18 to be par-mod, found '$ABI_ORDER'." >&2
    exit 1
fi

# ------------------------------------------------------------
# Exact group layout created by run_multi_molecule_MD.sh:
#
# 7  DTX modifier
# 8  DTX parent
# 9  Abi parent
# 10 Abi modifier
# 14 Pol
# 15 DTX-SI-C18
# 16 Abi-SI-C18
# 17 total Pro
# 18 Pol+Pro
# 19 PEG
# 20 PCL
# ------------------------------------------------------------
G_DTX_MOD=7
G_DTX_PAR=8
G_ABI_PAR=9
G_ABI_MOD=10

G_POL=14
G_DTX_PRO=15
G_ABI_PRO=16
G_PRO=17
G_POL_PRO=18
G_PEG=19
G_PCL=20

cat > "$OUT/group_mapping.log" <<EOF
DTX key = $DTX_KEY
DTX residue1/residue2/order = $DTX_RES1 / $DTX_RES2 / $DTX_ORDER
DTX Pro/Par/Mod = $G_DTX_PRO / $G_DTX_PAR / $G_DTX_MOD

Abi key = $ABI_KEY
Abi residue1/residue2/order = $ABI_RES1 / $ABI_RES2 / $ABI_ORDER
Abi Pro/Par/Mod = $G_ABI_PRO / $G_ABI_PAR / $G_ABI_MOD

Total Pro = $G_PRO
Pol = $G_POL
PEG = $G_PEG
PCL = $G_PCL
Pol+Pro = $G_POL_PRO
EOF

echo
echo "============================================================"
echo "Mixed DTX-SI-C18 + Abi-SI-C18 analysis"
echo "============================================================"
cat "$OUT/group_mapping.log"
echo

# ============================================================
# 1. PBC correction
# ============================================================
echo "[1/13] PBC correction"

printf "%s\n%s\n0\n" "$G_POL_PRO" "$G_POL_PRO" | \
"$GMX" trjconv \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$SYSTEM_DIR/MD.xtc" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/MD_cluster.xtc" \
    -pbc cluster \
    -center

printf "%s\n0\n" "$G_POL_PRO" | \
"$GMX" trjconv \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$OUT/MD_cluster.xtc" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/MD_noPBC.xtc" \
    -pbc whole \
    -center

rm -f "$OUT/MD_cluster.xtc"
TRAJ="$OUT/MD_noPBC.xtc"

# ============================================================
# 2. Energy QC
# ============================================================
echo "[2/13] Energy QC"

if [ -f "$SYSTEM_DIR/em.edr" ]; then
    printf "Potential\n0\n" | \
    "$GMX" energy -f "$SYSTEM_DIR/em.edr" -o "$OUT/potential.xvg"
fi

if [ -f "$SYSTEM_DIR/NVT.edr" ]; then
    printf "Temperature\n0\n" | \
    "$GMX" energy -f "$SYSTEM_DIR/NVT.edr" -o "$OUT/temperature.xvg"
fi

if [ -f "$SYSTEM_DIR/NPT.edr" ]; then
    printf "Pressure\n0\n" | \
    "$GMX" energy -f "$SYSTEM_DIR/NPT.edr" -o "$OUT/pressure.xvg"

    printf "Density\n0\n" | \
    "$GMX" energy -f "$SYSTEM_DIR/NPT.edr" -o "$OUT/density.xvg"
fi

# ============================================================
# 3. RMSD
# ============================================================
echo "[3/13] RMSD"

printf "%s\n%s\n" "$G_POL_PRO" "$G_POL_PRO" | \
"$GMX" rms -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/rmsd_polymer_drug.xvg" -tu ns

printf "%s\n%s\n" "$G_PRO" "$G_PRO" | \
"$GMX" rms -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/rmsd_drug.xvg" -tu ns

printf "%s\n%s\n" "$G_DTX_PRO" "$G_DTX_PRO" | \
"$GMX" rms -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/rmsd_dtxsic18.xvg" -tu ns

printf "%s\n%s\n" "$G_ABI_PRO" "$G_ABI_PRO" | \
"$GMX" rms -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/rmsd_abisic18.xvg" -tu ns

# ============================================================
# 4. Radius of gyration
# ============================================================
echo "[4/13] Radius of gyration"

printf "%s\n" "$G_POL_PRO" | \
"$GMX" gyrate -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/gyrate_polymer_drug.xvg"

printf "%s\n" "$G_POL" | \
"$GMX" gyrate -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/gyrate_polymer.xvg"

printf "%s\n" "$G_PRO" | \
"$GMX" gyrate -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/gyrate_total_drug.xvg"

printf "%s\n" "$G_DTX_PRO" | \
"$GMX" gyrate -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/gyrate_dtxsic18.xvg"

printf "%s\n" "$G_ABI_PRO" | \
"$GMX" gyrate -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/gyrate_abisic18.xvg"

printf "%s\n" "$G_PEG" | \
"$GMX" gyrate -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/gyrate_peg.xvg"

printf "%s\n" "$G_PCL" | \
"$GMX" gyrate -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/gyrate_pcl.xvg"

# ============================================================
# 5. COM distances: whole components
# ============================================================
echo "[5/13] Whole-component COM distances"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PRO plus com of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PRO plus com of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abi-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_PRO plus com of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-abi-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PRO plus com of group $G_ABI_PRO"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PRO plus com of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PRO plus com of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abi-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_PRO plus com of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PRO plus com of group $G_PEG"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PRO plus com of group $G_PEG"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abi-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_PRO plus com of group $G_PEG"

# ============================================================
# 6. COM distances: segment-resolved
# ============================================================
echo "[6/13] Segment-resolved COM distances"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxpar-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PAR plus com of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxlig-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_MOD plus com of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxpar-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PAR plus com of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxlig-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_MOD plus com of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxpar-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PAR plus com of group $G_PEG"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxlig-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_MOD plus com of group $G_PEG"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abipar-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_PAR plus com of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abilig-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_MOD plus com of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abipar-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_PAR plus com of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abilig-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_MOD plus com of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abipar-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_PAR plus com of group $G_PEG"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abilig-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_ABI_MOD plus com of group $G_PEG"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxpar-abipar-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PAR plus com of group $G_ABI_PAR"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxlig-abilig-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_MOD plus com of group $G_ABI_MOD"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxpar-abilig-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_PAR plus com of group $G_ABI_MOD"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtxlig-abipar-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_DTX_MOD plus com of group $G_ABI_PAR"

# ============================================================
# 7. COG distances used by the historical analysis
# ============================================================
echo "[7/13] COG distances"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-polymer-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_PRO plus cog of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-polymer-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_DTX_PRO plus cog of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abi-polymer-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_ABI_PRO plus cog of group $G_POL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-abi-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_DTX_PRO plus cog of group $G_ABI_PRO"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-pcl-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_DTX_PRO plus cog of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abi-pcl-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_ABI_PRO plus cog of group $G_PCL"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/dtx-peg-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_DTX_PRO plus cog of group $G_PEG"

"$GMX" distance -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/abi-peg-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_ABI_PRO plus cog of group $G_PEG"

# ============================================================
# 8. Hydrogen bonds
# ============================================================
echo "[8/13] Hydrogen bonds"

# Whole-component interactions
printf "%s\n%s\n" "$G_PRO" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-drug-polymer.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtx-polymer.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_PRO" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abi-polymer.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_ABI_PRO" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtx-abi.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtx-pcl.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_PRO" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abi-pcl.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtx-peg.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_PRO" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abi-peg.xvg" -dt 500

# DTX segment-carrier interactions
printf "%s\n%s\n" "$G_DTX_PAR" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxpar-polymer.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_MOD" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxlig-polymer.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_PAR" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxpar-pcl.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_MOD" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxlig-pcl.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_PAR" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxpar-peg.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_MOD" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxlig-peg.xvg" -dt 500

# Abi segment-carrier interactions
printf "%s\n%s\n" "$G_ABI_PAR" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abipar-polymer.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_MOD" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abilig-polymer.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_PAR" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abipar-pcl.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_MOD" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abilig-pcl.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_PAR" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abipar-peg.xvg" -dt 500

printf "%s\n%s\n" "$G_ABI_MOD" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-abilig-peg.xvg" -dt 500

# DTX-Abi segment-level interactions
printf "%s\n%s\n" "$G_DTX_PAR" "$G_ABI_PAR" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxpar-abipar.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_MOD" "$G_ABI_MOD" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxlig-abilig.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_PAR" "$G_ABI_MOD" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxpar-abilig.xvg" -dt 500

printf "%s\n%s\n" "$G_DTX_MOD" "$G_ABI_PAR" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dtxlig-abipar.xvg" -dt 500

# ============================================================
# 9. Minimum distances and contact numbers
# ============================================================
echo "[9/13] Minimum distances and contacts"

printf "%s\n%s\n" "$G_DTX_PRO" "$G_ABI_PRO" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/dtx-abi-mindist.xvg" -on "$OUT/dtx-abi-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_POL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/dtx-polymer-mindist.xvg" -on "$OUT/dtx-polymer-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_ABI_PRO" "$G_POL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/abi-polymer-mindist.xvg" -on "$OUT/abi-polymer-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_PCL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/dtx-pcl-mindist.xvg" -on "$OUT/dtx-pcl-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_ABI_PRO" "$G_PCL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/abi-pcl-mindist.xvg" -on "$OUT/abi-pcl-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_PEG" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/dtx-peg-mindist.xvg" -on "$OUT/dtx-peg-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_ABI_PRO" "$G_PEG" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/abi-peg-mindist.xvg" -on "$OUT/abi-peg-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_DTX_MOD" "$G_PCL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/dtxlig-pcl-mindist.xvg" -on "$OUT/dtxlig-pcl-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_DTX_PAR" "$G_PCL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/dtxpar-pcl-mindist.xvg" -on "$OUT/dtxpar-pcl-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_ABI_MOD" "$G_PCL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/abilig-pcl-mindist.xvg" -on "$OUT/abilig-pcl-contacts.xvg" \
    -d 0.6 -dt 500

printf "%s\n%s\n" "$G_ABI_PAR" "$G_PCL" | \
"$GMX" mindist -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -od "$OUT/abipar-pcl-mindist.xvg" -on "$OUT/abipar-pcl-contacts.xvg" \
    -d 0.6 -dt 500

# ============================================================
# 10. RDF
# Molecular centres of mass are used for both reference and selected groups.
# ============================================================
echo "[10/13] RDF"

"$GMX" rdf -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -ref "group $G_DTX_PRO" -sel "group $G_PCL" \
    -selrpos mol_com -seltype mol_com -pbc yes -bin 0.002 \
    -o "$OUT/rdf-dtx-pcl.xvg" -cn "$OUT/rdf-dtx-pcl-cn.xvg"

"$GMX" rdf -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -ref "group $G_ABI_PRO" -sel "group $G_PCL" \
    -selrpos mol_com -seltype mol_com -pbc yes -bin 0.002 \
    -o "$OUT/rdf-abi-pcl.xvg" -cn "$OUT/rdf-abi-pcl-cn.xvg"

"$GMX" rdf -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -ref "group $G_DTX_PRO" -sel "group $G_PEG" \
    -selrpos mol_com -seltype mol_com -pbc yes -bin 0.002 \
    -o "$OUT/rdf-dtx-peg.xvg" -cn "$OUT/rdf-dtx-peg-cn.xvg"

"$GMX" rdf -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -ref "group $G_ABI_PRO" -sel "group $G_PEG" \
    -selrpos mol_com -seltype mol_com -pbc yes -bin 0.002 \
    -o "$OUT/rdf-abi-peg.xvg" -cn "$OUT/rdf-abi-peg-cn.xvg"

"$GMX" rdf -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -ref "group $G_DTX_PRO" -sel "group $G_ABI_PRO" \
    -selrpos mol_com -seltype mol_com -pbc yes -bin 0.002 \
    -o "$OUT/rdf-dtx-abi.xvg" -cn "$OUT/rdf-dtx-abi-cn.xvg"

"$GMX" rdf -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -ref "group $G_PRO" -sel "group $G_PCL" \
    -selrpos mol_com -seltype mol_com -pbc yes -bin 0.002 \
    -o "$OUT/rdf-drug-pcl.xvg" -cn "$OUT/rdf-drug-pcl-cn.xvg"

"$GMX" rdf -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -ref "group $G_PRO" -sel "group $G_PEG" \
    -selrpos mol_com -seltype mol_com -pbc yes -bin 0.002 \
    -o "$OUT/rdf-drug-peg.xvg" -cn "$OUT/rdf-drug-peg-cn.xvg"

# ============================================================
# 11. SASA
# Contact area can be calculated as:
# (SASA_A + SASA_B - SASA_A+B) / 2
# ============================================================
echo "[11/13] SASA"

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_POL" -o "$OUT/sasa-polymer.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PRO" -o "$OUT/sasa-drug.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_DTX_PRO" -o "$OUT/sasa-dtx.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_ABI_PRO" -o "$OUT/sasa-abi.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PCL" -o "$OUT/sasa-pcl.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PEG" -o "$OUT/sasa-peg.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_DTX_PRO or group $G_POL" \
    -o "$OUT/sasa-dtx-polymer.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_ABI_PRO or group $G_POL" \
    -o "$OUT/sasa-abi-polymer.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PRO or group $G_POL" \
    -o "$OUT/sasa-drug-polymer.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_DTX_PRO or group $G_ABI_PRO" \
    -o "$OUT/sasa-dtx-abi.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_DTX_PRO or group $G_PCL" \
    -o "$OUT/sasa-dtx-pcl.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_ABI_PRO or group $G_PCL" \
    -o "$OUT/sasa-abi-pcl.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_DTX_PRO or group $G_PEG" \
    -o "$OUT/sasa-dtx-peg.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_ABI_PRO or group $G_PEG" \
    -o "$OUT/sasa-abi-peg.xvg" -tu ns -dt 0.5

# ============================================================
# 12. Cluster analysis
# ============================================================
echo "[12/13] Cluster analysis"

printf "%s\n%s\n" "$G_PRO" "$G_PRO" | \
"$GMX" cluster \
    -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -cutoff 2 \
    -o "$OUT/cluster_drug.xpm" \
    -g "$OUT/cluster_drug.log" \
    -sz "$OUT/size_drug.xvg" \
    -ntr "$OUT/ntr_drug.xvg" \
    -clid "$OUT/nclid_drug.xvg" \
    -cl "$OUT/clusters_drug.pdb" \
    -tu ps -dt 500

printf "%s\n%s\n" "$G_DTX_PRO" "$G_DTX_PRO" | \
"$GMX" cluster \
    -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -cutoff 2 \
    -o "$OUT/cluster_dtx.xpm" \
    -g "$OUT/cluster_dtx.log" \
    -sz "$OUT/size_dtx.xvg" \
    -ntr "$OUT/ntr_dtx.xvg" \
    -clid "$OUT/nclid_dtx.xvg" \
    -cl "$OUT/clusters_dtx.pdb" \
    -tu ps -dt 500

printf "%s\n%s\n" "$G_ABI_PRO" "$G_ABI_PRO" | \
"$GMX" cluster \
    -f "$TRAJ" -s "$SYSTEM_DIR/MD.tpr" -n "$SYSTEM_DIR/index.ndx" \
    -cutoff 2 \
    -o "$OUT/cluster_abi.xpm" \
    -g "$OUT/cluster_abi.log" \
    -sz "$OUT/size_abi.xvg" \
    -ntr "$OUT/ntr_abi.xvg" \
    -clid "$OUT/nclid_abi.xvg" \
    -cl "$OUT/clusters_abi.pdb" \
    -tu ps -dt 500

# ============================================================
# 13. Representative trajectory
# ============================================================
echo "[13/13] Representative trajectory"

printf "%s\n" "$G_POL_PRO" | \
"$GMX" trjconv \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$TRAJ" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/movie-polymer-drug-noPBC.pdb" \
    -dt 500

echo
echo "DONE"
echo "Output: $OUT"
