#!/usr/bin/env bash
set -euo pipefail

# Self-association MM-PBSA.
#
# Usage:
#   bash mmpbsa_analysis.sh PATH_TO_SYSTEM MOL_KEY
#
# Run prepare_mmpbsa_index.sh first.

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 PATH_TO_SYSTEM MOL_KEY"
    exit 1
fi

NPROC="${NPROC:-6}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$POST_ROOT/.." && pwd)"

SEGMENT_TABLE="$POST_ROOT/segment_definition.tsv"
MD_ROOT="$REPO_ROOT/MD_simulation/prodrug_self_association_MD"
OUTPUT_ROOT="$POST_ROOT/analysis_output/prodrug_self_association_MD"
MMPBSA_INPUT="$REPO_ROOT/MD_simulation/MDP_files/mmpbsa.in"

SYSTEM_DIR="$(cd "$1" && pwd)"
MOL_KEY="$2"

case "$SYSTEM_DIR/" in
    "$MD_ROOT/"*) RELATIVE_SYSTEM="${SYSTEM_DIR#"$MD_ROOT/"}" ;;
    *) echo "ERROR: system is not inside $MD_ROOT" >&2; exit 1 ;;
esac

OUT="$OUTPUT_ROOT/$RELATIVE_SYSTEM"
TRAJ="$OUT/MD_noPBC.xtc"
MMPBSA_INDEX="$OUT/index_mmpbsa.ndx"
MMPBSA_OUT="$OUT/mmpbsa"
mkdir -p "$MMPBSA_OUT"

for FILE in MD.tpr topol.top; do
    [ -f "$SYSTEM_DIR/$FILE" ] || { echo "ERROR: missing $SYSTEM_DIR/$FILE" >&2; exit 1; }
done

[ -f "$TRAJ" ] || { echo "ERROR: run prepare_mmpbsa_index.sh first." >&2; exit 1; }
[ -f "$MMPBSA_INDEX" ] || { echo "ERROR: run prepare_mmpbsa_index.sh first." >&2; exit 1; }
[ -f "$MMPBSA_INPUT" ] || { echo "ERROR: missing $MMPBSA_INPUT" >&2; exit 1; }


# Read one manually verified row from segment_definition.tsv.
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


echo
echo "============================================================"
echo "Self-association MM-PBSA: $RELATIVE_SYSTEM"
echo "mol = $MOL_KEY"
echo "order = $ORDER"
echo "============================================================"

cd "$MMPBSA_OUT"

if [ "$ORDER" = "par" ]; then
    # index_mmpbsa.ndx:
    # 7 = Pro1
    # 8 = Pro2
    echo "[1/1] Pro-Pro"

    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$MMPBSA_INDEX" \
        -cg 7 8 \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-Pro.dat \
        -eo Pro-Pro.csv \
        -nogui
else
    # The generated prodrug index always has:
    # 8/9   = residue1 molecule1/molecule2
    # 10/11 = residue2 molecule1/molecule2
    # 12/13 = Pro1/Pro2

    if [ "$ORDER" = "mod-par" ]; then
        G_MOD1=8
        G_MOD2=9
        G_PAR1=10
        G_PAR2=11
    else
        G_PAR1=8
        G_PAR2=9
        G_MOD1=10
        G_MOD2=11
    fi

    G_PRO1=12
    G_PRO2=13

    echo "[1/5] Pro-Pro"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$MMPBSA_INDEX" \
        -cg "$G_PRO1" "$G_PRO2" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Pro-Pro.dat \
        -eo Pro-Pro.csv \
        -nogui

    echo "[2/5] Par-Par"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$MMPBSA_INDEX" \
        -cg "$G_PAR1" "$G_PAR2" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Par-Par.dat \
        -eo Par-Par.csv \
        -nogui

    echo "[3/5] Mod-Mod"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$MMPBSA_INDEX" \
        -cg "$G_MOD1" "$G_MOD2" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Mod-Mod.dat \
        -eo Mod-Mod.csv \
        -nogui

    echo "[4/5] Par1-Mod2"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$MMPBSA_INDEX" \
        -cg "$G_PAR1" "$G_MOD2" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Par1-Mod2.dat \
        -eo Par1-Mod2.csv \
        -nogui

    echo "[5/5] Par2-Mod1"
    mpirun -np "$NPROC" gmx_MMPBSA \
        -O -i "$MMPBSA_INPUT" \
        -cs "$SYSTEM_DIR/MD.tpr" \
        -ci "$MMPBSA_INDEX" \
        -cg "$G_PAR2" "$G_MOD1" \
        -ct "$TRAJ" \
        -cp "$SYSTEM_DIR/topol.top" \
        -o Par2-Mod1.dat \
        -eo Par2-Mod1.csv \
        -nogui

    # SCRAP Par-Mod descriptor: arithmetic mean of the two reciprocal
    # cross-molecule interactions, Par1-Mod2 and Par2-Mod1.
    extract_delta_total() {
        local file="$1"
        awk '''$1 == "DELTA" && $2 == "TOTAL" { value=$3 }
             END { if (value == "") exit 1; print value }''' "$file"
    }

    PAR1_MOD2="$(extract_delta_total Par1-Mod2.dat)" || {
        echo "ERROR: could not extract DELTA TOTAL from Par1-Mod2.dat" >&2
        exit 1
    }
    PAR2_MOD1="$(extract_delta_total Par2-Mod1.dat)" || {
        echo "ERROR: could not extract DELTA TOTAL from Par2-Mod1.dat" >&2
        exit 1
    }
    PAR_MOD_MEAN="$(awk -v a="$PAR1_MOD2" -v b="$PAR2_MOD1" 'BEGIN { printf "%.6f", (a+b)/2 }')"

    {
        printf "term\tDeltaG_kcal_mol-1\n"
        printf "Par1-Mod2\t%s\n" "$PAR1_MOD2"
        printf "Par2-Mod1\t%s\n" "$PAR2_MOD1"
        printf "Par-Mod\t%s\n" "$PAR_MOD_MEAN"
    } > Par-Mod_summary.tsv

    echo "Par-Mod = mean(Par1-Mod2, Par2-Mod1) = $PAR_MOD_MEAN kcal/mol"
fi

echo
echo "DONE"
echo "Output: $MMPBSA_OUT"
