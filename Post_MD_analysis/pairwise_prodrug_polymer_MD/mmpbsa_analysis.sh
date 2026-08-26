#!/usr/bin/env bash
set -euo pipefail

# Pairwise carrier-facing MM-PBSA.
#
# Usage:
#   bash mmpbsa_analysis.sh PATH_TO_SYSTEM MOL_KEY
#
# Run gromacs_analysis.sh first.
#
# IMPORTANT: -cg follows the original deposited workflow:
# carrier group first, drug/prodrug segment second.
# Example: Par-PCL uses -cg "$G_PCL" "$G_PAR".

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 PATH_TO_SYSTEM MOL_KEY"
    exit 1
fi

NPROC="${NPROC:-6}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$POST_ROOT/.." && pwd)"

SEGMENT_TABLE="$POST_ROOT/segment_defination.tsv"
MD_ROOT="$REPO_ROOT/MD_simulation/pairwise_prodrug_polymer_MD"
OUTPUT_ROOT="$POST_ROOT/analysis_output/pairwise_prodrug_polymer_MD"
MMPBSA_INPUT="$REPO_ROOT/MD_simulation/MDP_files/mmpbsa.in"

SYSTEM_DIR="$(cd "$1" && pwd)"
MOL_KEY="$2"

case "$SYSTEM_DIR/" in
    "$MD_ROOT/"*) RELATIVE_SYSTEM="${SYSTEM_DIR#"$MD_ROOT/"}" ;;
    *) echo "ERROR: system is not inside $MD_ROOT" >&2; exit 1 ;;
esac

OUT="$OUTPUT_ROOT/$RELATIVE_SYSTEM"
TRAJ="$OUT/MD_noPBC.xtc"
MMPBSA_OUT="$OUT/mmpbsa"
mkdir -p "$MMPBSA_OUT"

for FILE in MD.tpr index.ndx topol.top; do
    [ -f "$SYSTEM_DIR/$FILE" ] || { echo "ERROR: missing $SYSTEM_DIR/$FILE" >&2; exit 1; }
done

[ -f "$TRAJ" ] || { echo "ERROR: run gromacs_analysis.sh first." >&2; exit 1; }
[ -f "$MMPBSA_INPUT" ] || { echo "ERROR: missing $MMPBSA_INPUT" >&2; exit 1; }


# Read one manually verified row from segment_defination.tsv.
# The second command-line argument is the value in the `mol` column.
ROW="$(
    awk -F '\t' -v key="$MOL_KEY" '
        NR > 1 && $1 == key {
            print $2 "|" $3 "|" $4
        }
    ' "$SEGMENT_TABLE"
)"

if [ -z "$ROW" ]; then
    echo "ERROR: '$MOL_KEY' was not found in $SEGMENT_TABLE" >&2
    exit 1
fi

if [ "$(printf '%s\n' "$ROW" | wc -l)" -ne 1 ]; then
    echo "ERROR: '$MOL_KEY' occurs more than once in $SEGMENT_TABLE" >&2
    exit 1
fi

IFS='|' read -r RESIDUE1 RESIDUE2 ORDER <<< "$ROW"

case "$ORDER" in
    par)
        PAR_RESIDUE="$RESIDUE1"
        MOD_RESIDUE=""
        ;;
    par-mod)
        PAR_RESIDUE="$RESIDUE1"
        MOD_RESIDUE="$RESIDUE2"
        ;;
    mod-par)
        MOD_RESIDUE="$RESIDUE1"
        PAR_RESIDUE="$RESIDUE2"
        ;;
    *)
        echo "ERROR: unsupported order '$ORDER' for $MOL_KEY" >&2
        exit 1
        ;;
esac

# Print the name of one numerical group from index.ndx.
index_group_name() {
    local wanted="$1"

    awk -v wanted="$wanted" '
        /^[[:space:]]*\[/ {
            name=$0
            sub(/^[[:space:]]*\[[[:space:]]*/, "", name)
            sub(/[[:space:]]*\][[:space:]]*$/, "", name)

            if (group_number == wanted) {
                print name
                exit
            }

            group_number++
        }
    ' group_number=0 "$SYSTEM_DIR/index.ndx"
}

check_group() {
    local number="$1"
    local expected="$2"
    local actual

    actual="$(index_group_name "$number")"

    if [ "$actual" != "$expected" ]; then
        echo "ERROR: index group $number is '$actual'; expected '$expected'." >&2
        exit 1
    fi
}


# Exact group layout from the MD batch script.
if [ "$ORDER" = "par" ]; then
    G_PRO=7
    G_PAR=7
    G_MOD=""
    G_POL=11
    G_PEG=13
    G_PCL=14

    check_group 7  "$RESIDUE1"
    check_group 11 "polymer"
    check_group 13 "peg"
    check_group 14 "pcl"
else
    G_PRO=13
    G_POL=12
    G_PEG=15
    G_PCL=16

    check_group 7  "$RESIDUE1"
    check_group 8  "$RESIDUE2"
    check_group 12 "polymer"
    check_group 13 "drug"
    check_group 15 "peg"
    check_group 16 "pcl"

    if [ "$ORDER" = "mod-par" ]; then
        G_MOD=7
        G_PAR=8
    else
        G_PAR=7
        G_MOD=8
    fi
fi

echo
echo "============================================================"
echo "MM-PBSA: $RELATIVE_SYSTEM"
echo "mol = $MOL_KEY"
echo "order = $ORDER"
echo "Pro=$G_PRO  Par=$G_PAR  Mod=${G_MOD:-NA}"
echo "Pol=$G_POL  PEG=$G_PEG  PCL=$G_PCL"
echo "============================================================"

cd "$MMPBSA_OUT"

if [ "$ORDER" = "par" ]; then
    echo "[1/3] Pro-Pol"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_POL" "$G_PRO" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-Pol.dat \
        -eo Pro-Pol.csv \
        -nogui

    echo "[2/3] Pro-PCL"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PCL" "$G_PRO" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-PCL.dat \
        -eo Pro-PCL.csv \
        -nogui

    echo "[3/3] Pro-PEG"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PEG" "$G_PRO" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-PEG.dat \
        -eo Pro-PEG.csv \
        -nogui
else
    echo "[1/9] Pro-Pol"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_POL" "$G_PRO" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-Pol.dat \
        -eo Pro-Pol.csv \
        -nogui

    echo "[2/9] Par-Pol"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_POL" "$G_PAR" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Par-Pol.dat \
        -eo Par-Pol.csv \
        -nogui

    echo "[3/9] Mod-Pol"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_POL" "$G_MOD" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Mod-Pol.dat \
        -eo Mod-Pol.csv \
        -nogui

    echo "[4/9] Pro-PCL"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PCL" "$G_PRO" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-PCL.dat \
        -eo Pro-PCL.csv \
        -nogui

    echo "[5/9] Par-PCL"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PCL" "$G_PAR" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Par-PCL.dat \
        -eo Par-PCL.csv \
        -nogui

    echo "[6/9] Mod-PCL"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PCL" "$G_MOD" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Mod-PCL.dat \
        -eo Mod-PCL.csv \
        -nogui

    echo "[7/9] Pro-PEG"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PEG" "$G_PRO" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-PEG.dat \
        -eo Pro-PEG.csv \
        -nogui

    echo "[8/9] Par-PEG"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PEG" "$G_PAR" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Par-PEG.dat \
        -eo Par-PEG.csv \
        -nogui

    echo "[9/9] Mod-PEG"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$SYSTEM_DIR/index.ndx" \
        -cg "$G_PEG" "$G_MOD" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Mod-PEG.dat \
        -eo Mod-PEG.csv \
        -nogui
fi

echo
echo "DONE"
echo "Output: $MMPBSA_OUT"
