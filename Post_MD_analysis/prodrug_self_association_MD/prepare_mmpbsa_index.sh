#!/usr/bin/env bash
set -euo pipefail

# Prepare the molecule-specific MM-PBSA index for a self-association system.
# No ITP file is read.
#
# Usage:
#   bash prepare_mmpbsa_index.sh PATH_TO_SYSTEM MOL_KEY

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 PATH_TO_SYSTEM MOL_KEY"
    exit 1
fi

GMX="${GMX:-gmx}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$POST_ROOT/.." && pwd)"

SEGMENT_TABLE="$POST_ROOT/segment_defination.tsv"
MD_ROOT="$REPO_ROOT/MD_simulation/prodrug_self_association_MD"
OUTPUT_ROOT="$POST_ROOT/analysis_output/prodrug_self_association_MD"

SYSTEM_DIR="$(cd "$1" && pwd)"
MOL_KEY="$2"

case "$SYSTEM_DIR/" in
    "$MD_ROOT/"*) RELATIVE_SYSTEM="${SYSTEM_DIR#"$MD_ROOT/"}" ;;
    *) echo "ERROR: system is not inside $MD_ROOT" >&2; exit 1 ;;
esac

OUT="$OUTPUT_ROOT/$RELATIVE_SYSTEM"
mkdir -p "$OUT"

for FILE in MD.tpr MD.xtc MD.gro index.ndx; do
    [ -f "$SYSTEM_DIR/$FILE" ] || { echo "ERROR: missing $SYSTEM_DIR/$FILE" >&2; exit 1; }
done


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


echo
echo "============================================================"
echo "Self-association preparation: $RELATIVE_SYSTEM"
echo "mol = $MOL_KEY"
echo "residue1 = $RESIDUE1"
echo "residue2 = ${RESIDUE2:-NA}"
echo "order = $ORDER"
echo "============================================================"

# -------------------------------------------------------------------------
# 1. Validate the original simulation index.
# -------------------------------------------------------------------------
if [ "$ORDER" = "par" ]; then
    check_group 2 "$RESIDUE1"
    check_group 6 "drug"
    G_DRUG=6
else
    check_group 2 "$RESIDUE1"
    check_group 3 "$RESIDUE2"
    check_group 7 "drug"
    G_DRUG=7
fi

# -------------------------------------------------------------------------
# 2. PBC correction of the two-drug system.
# -------------------------------------------------------------------------
echo "[1] PBC correction"

printf "%s\n%s\n0\n" "$G_DRUG" "$G_DRUG" | \
"$GMX" trjconv \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$SYSTEM_DIR/MD.xtc" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/MD_cluster.xtc" \
    -pbc cluster \
    -center

printf "%s\n0\n" "$G_DRUG" | \
"$GMX" trjconv \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$OUT/MD_cluster.xtc" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/MD_noPBC.xtc" \
    -pbc whole \
    -center

rm -f "$OUT/MD_cluster.xtc"

# -------------------------------------------------------------------------
# 3. Create molecule-specific index groups with gmx make_ndx.
#
# The original self-association prodrug index ends at group 7.
# Therefore:
#
# splitres 2 -> new groups 8 and 9
# splitres 3 -> new groups 10 and 11
# 8|10       -> new group 12 = Pro1
# 9|11       -> new group 13 = Pro2
#
# The TSV determines whether groups 8/9 or 10/11 are Par/Mod.
# -------------------------------------------------------------------------
echo "[2] Create molecule-specific MM-PBSA index"

if [ "$ORDER" = "par" ]; then
    # Original parent-drug index ends at group 6.
    # splitres 2 creates groups 7 and 8, one for each parent-drug molecule.
    printf "splitres 2\nname 7 Pro1\nname 8 Pro2\nq\n" | \
    "$GMX" make_ndx \
        -f "$SYSTEM_DIR/MD.gro" \
        -n "$SYSTEM_DIR/index.ndx" \
        -o "$OUT/index_mmpbsa.ndx"

    cat > "$OUT/group_mapping.log" <<EOF
mol = $MOL_KEY
order = par

Original simulation:
  group 2 = $RESIDUE1
  group 6 = drug

Generated MM-PBSA index:
  group 7 = Pro1
  group 8 = Pro2
EOF
else
    printf "splitres 2\nsplitres 3\n8|10\nname 12 Pro1\n9|11\nname 13 Pro2\nq\n" | \
    "$GMX" make_ndx \
        -f "$SYSTEM_DIR/MD.gro" \
        -n "$SYSTEM_DIR/index.ndx" \
        -o "$OUT/index_mmpbsa.ndx"

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

    cat > "$OUT/group_mapping.log" <<EOF
mol = $MOL_KEY
residue1 = $RESIDUE1
residue2 = $RESIDUE2
order = $ORDER

Original simulation:
  group 2 = residue1
  group 3 = residue2
  group 7 = drug

Generated by splitres:
  group 8  = residue1_mol1
  group 9  = residue1_mol2
  group 10 = residue2_mol1
  group 11 = residue2_mol2

SCRAP groups:
  Par1 = $G_PAR1
  Par2 = $G_PAR2
  Mod1 = $G_MOD1
  Mod2 = $G_MOD2
  Pro1 = 12
  Pro2 = 13
EOF
fi

echo
cat "$OUT/group_mapping.log"
echo
echo "DONE"
echo "Output: $OUT"
