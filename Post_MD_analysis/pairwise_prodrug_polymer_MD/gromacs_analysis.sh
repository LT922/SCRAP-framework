#!/usr/bin/env bash
set -euo pipefail

# Pairwise drug/prodrug-polymer GROMACS post-processing.
#
# Usage:
#   bash gromacs_analysis.sh PATH_TO_SYSTEM MOL_KEY
#
# Example:
#   bash gromacs_analysis.sh \
#     ../../MD_simulation/pairwise_prodrug_polymer_MD/training/P2-C18 \
#     ptx2c18

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 PATH_TO_SYSTEM MOL_KEY"
    exit 1
fi

GMX="${GMX:-gmx}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$POST_ROOT/.." && pwd)"

SEGMENT_TABLE="$POST_ROOT/segment_definition.tsv"
MD_ROOT="$REPO_ROOT/MD_simulation/pairwise_prodrug_polymer_MD"
OUTPUT_ROOT="$POST_ROOT/analysis_output/pairwise_prodrug_polymer_MD"

SYSTEM_DIR="$(cd "$1" && pwd)"
MOL_KEY="$2"

case "$SYSTEM_DIR/" in
    "$MD_ROOT/"*)
        RELATIVE_SYSTEM="${SYSTEM_DIR#"$MD_ROOT/"}"
        ;;
    *)
        echo "ERROR: system is not inside $MD_ROOT" >&2
        exit 1
        ;;
esac

OUT="$OUTPUT_ROOT/$RELATIVE_SYSTEM"
mkdir -p "$OUT"

for FILE in MD.tpr MD.xtc index.ndx; do
    if [ ! -f "$SYSTEM_DIR/$FILE" ]; then
        echo "ERROR: missing $SYSTEM_DIR/$FILE" >&2
        exit 1
    fi
done


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


# -------------------------------------------------------------------------
# 1. Convert the known MD index layout into Pro/Par/Mod/Pol/PEG/PCL.
# -------------------------------------------------------------------------
if [ "$ORDER" = "par" ]; then
    # Parent-drug index:
    # 7 parent drug; 11 polymer; 12 polymer&drug; 13 peg; 14 pcl
    G_PRO=7
    G_PAR=7
    G_MOD=""
    G_POL=11
    G_POL_PRO=12
    G_PEG=13
    G_PCL=14

    check_group 7  "$RESIDUE1"
    check_group 11 "polymer"
    check_group 12 "polymer&drug"
    check_group 13 "peg"
    check_group 14 "pcl"
else
    # Prodrug index:
    # 7 residue1; 8 residue2; 12 polymer; 13 drug;
    # 14 polymer&drug; 15 peg; 16 pcl
    G_PRO=13
    G_POL=12
    G_POL_PRO=14
    G_PEG=15
    G_PCL=16

    check_group 7  "$RESIDUE1"
    check_group 8  "$RESIDUE2"
    check_group 12 "polymer"
    check_group 13 "drug"
    check_group 14 "polymer&drug"
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

cat > "$OUT/group_mapping.log" <<EOF
mol = $MOL_KEY
residue1 = $RESIDUE1
residue2 = ${RESIDUE2:-NA}
order = $ORDER

Pro = $G_PRO
Par = $G_PAR
Mod = ${G_MOD:-NA}
Pol = $G_POL
PEG = $G_PEG
PCL = $G_PCL
Pol+Pro = $G_POL_PRO
EOF

echo
echo "============================================================"
echo "Pairwise analysis: $RELATIVE_SYSTEM"
echo "============================================================"
cat "$OUT/group_mapping.log"
echo

# -------------------------------------------------------------------------
# 2. PBC correction.
# The Pol+Pro assembly is centred. Group 0 (System) is retained in the
# output trajectory, so all atom numbers remain compatible with MD.tpr
# and index.ndx.
# -------------------------------------------------------------------------
echo "[1] PBC correction"

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

# -------------------------------------------------------------------------
# 3. RMSD and radius of gyration.
# -------------------------------------------------------------------------
echo "[2] RMSD and Rg"

printf "%s\n%s\n" "$G_POL_PRO" "$G_POL_PRO" | \
"$GMX" rms \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$TRAJ" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/rmsd.xvg" \
    -tu ns

printf "%s\n" "$G_POL_PRO" | \
"$GMX" gyrate \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$TRAJ" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/allgyrate.xvg"

printf "%s\n" "$G_PRO" | \
"$GMX" gyrate \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$TRAJ" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/druggyrate.xvg"

printf "%s\n" "$G_POL" | \
"$GMX" gyrate \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$TRAJ" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/polymergyrate.xvg"

# -------------------------------------------------------------------------
# 4. COM distances.
# -------------------------------------------------------------------------
echo "[3] COM distances"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PRO plus com of group $G_POL"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/par-polymer-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PAR plus com of group $G_POL"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PRO plus com of group $G_PCL"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/par-pcl-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PAR plus com of group $G_PCL"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PRO plus com of group $G_PEG"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/par-peg-dist.xvg" -tu ns -dt 0.5 \
    -select "com of group $G_PAR plus com of group $G_PEG"

if [ -n "$G_MOD" ]; then
    "$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -oall "$OUT/lig-polymer-dist.xvg" -tu ns -dt 0.5 \
        -select "com of group $G_MOD plus com of group $G_POL"

    "$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -oall "$OUT/lig-pcl-dist.xvg" -tu ns -dt 0.5 \
        -select "com of group $G_MOD plus com of group $G_PCL"

    "$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -oall "$OUT/lig-peg-dist.xvg" -tu ns -dt 0.5 \
        -select "com of group $G_MOD plus com of group $G_PEG"
fi

# -------------------------------------------------------------------------
# 5. COG distances.
# -------------------------------------------------------------------------
echo "[4] COG distances"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/drug-polymer-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_PRO plus cog of group $G_POL"

"$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -oall "$OUT/par-polymer-cogdist.xvg" -tu ns -dt 0.5 \
    -select "cog of group $G_PAR plus cog of group $G_POL"

if [ -n "$G_MOD" ]; then
    "$GMX" distance -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -oall "$OUT/lig-polymer-cogdist.xvg" -tu ns -dt 0.5 \
        -select "cog of group $G_MOD plus cog of group $G_POL"
fi

# -------------------------------------------------------------------------
# 6. Hydrogen bonds.
#
# This is intentionally written in the same explicit form as the original
# workflow. Example for Par-PEG:
#
#   printf "$G_PAR \n $G_PEG \n" | gmx hbond ...
# -------------------------------------------------------------------------
echo "[5] Hydrogen bonds"

printf "%s\n%s\n" "$G_PRO" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dru-po.xvg" -tu ns -dt 0.5

printf "%s\n%s\n" "$G_PAR" "$G_POL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-par-po.xvg" -tu ns -dt 0.5

printf "%s\n%s\n" "$G_PRO" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dru-pcl.xvg" -tu ns -dt 0.5

printf "%s\n%s\n" "$G_PAR" "$G_PCL" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-par-pcl.xvg" -tu ns -dt 0.5

printf "%s\n%s\n" "$G_PRO" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-dru-peg.xvg" -tu ns -dt 0.5

printf "%s\n%s\n" "$G_PAR" "$G_PEG" | \
"$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -num "$OUT/hbnum-par-peg.xvg" -tu ns -dt 0.5

if [ -n "$G_MOD" ]; then
    printf "%s\n%s\n" "$G_MOD" "$G_POL" | \
    "$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -num "$OUT/hbnum-lig-po.xvg" -tu ns -dt 0.5

    printf "%s\n%s\n" "$G_MOD" "$G_PCL" | \
    "$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -num "$OUT/hbnum-lig-pcl.xvg" -tu ns -dt 0.5

    printf "%s\n%s\n" "$G_MOD" "$G_PEG" | \
    "$GMX" hbond -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -num "$OUT/hbnum-lig-peg.xvg" -tu ns -dt 0.5
fi

# -------------------------------------------------------------------------
# 7. SASA.
# Contact area between Pro and Pol can be calculated as:
#
#   (SASA_Pro + SASA_Pol - SASA_Pro+Pol) / 2
# -------------------------------------------------------------------------
echo "[6] SASA"

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PRO" -o "$OUT/area_dru.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PAR" -o "$OUT/area_par.xvg" -tu ns -dt 0.5

if [ -n "$G_MOD" ]; then
    "$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
        -surface "group $G_MOD" -o "$OUT/area_lig.xvg" -tu ns -dt 0.5
fi

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_POL" -o "$OUT/area_pol.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PCL" -o "$OUT/area_core.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_PEG" -o "$OUT/area_peg.xvg" -tu ns -dt 0.5

"$GMX" sasa -s "$SYSTEM_DIR/MD.tpr" -f "$TRAJ" -n "$SYSTEM_DIR/index.ndx" \
    -surface "group $G_POL_PRO" -o "$OUT/area_whole.xvg" -tu ns -dt 0.5

# -------------------------------------------------------------------------
# 8. Cluster analysis and representative trajectory.
# -------------------------------------------------------------------------
echo "[7] Cluster analysis"

printf "%s\n%s\n" "$G_POL_PRO" "$G_POL_PRO" | \
"$GMX" cluster \
    -f "$TRAJ" \
    -s "$SYSTEM_DIR/MD.tpr" \
    -n "$SYSTEM_DIR/index.ndx" \
    -cutoff 2 \
    -o "$OUT/cluster.xpm" \
    -g "$OUT/cluster.log" \
    -sz "$OUT/size.xvg" \
    -ntr "$OUT/ntr.xvg" \
    -clid "$OUT/nclid.xvg" \
    -cl "$OUT/clusters.pdb" \
    -tu ps \
    -dt 500

echo "[8] Representative trajectory"

printf "%s\n" "$G_POL_PRO" | \
"$GMX" trjconv \
    -s "$SYSTEM_DIR/MD.tpr" \
    -f "$TRAJ" \
    -n "$SYSTEM_DIR/index.ndx" \
    -o "$OUT/movie-noPBC.pdb" \
    -dt 500

echo
echo "DONE"
echo "Output: $OUT"
