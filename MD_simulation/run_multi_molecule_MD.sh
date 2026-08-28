#!/bin/bash

# ============================================================
# Batch MD: multi_molecule_MD
# Restart order:
# MD.gro -> skip
# MD.tpr + MD.cpt -> continue MD directly
# NPT.gro + NPT.cpt -> MD
# NVT.gro + NVT.cpt -> NPT -> MD
# em.gro -> NVT -> NPT -> MD
# sol.gro -> index -> EM -> NVT -> NPT -> MD
# polymerdrug.gro -> solvate -> ...
# polymerdrug.pdb -> editconf -> ...
# otherwise -> Packmol -> ...
# ============================================================

GMX="${GMX:-gmx}"
NTOMP="${NTOMP:-8}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/multi_molecule_MD"

if [ ! -d "$ROOT" ]; then
    echo "ERROR: multi_molecule_MD not found."
    exit 1
fi

: > "$ROOT/MDdone.log"
: > "$ROOT/MDfailed.log"

for folder in "$ROOT"/*
do
    [ -d "$folder" ] || continue
    [ -f "$folder/packmol.inp" ] || continue
    [ -f "$folder/topol.top" ] || continue

    name="$(basename "$folder")"

        echo
        echo "============================================================"
        echo "Processing: $name"
        echo "============================================================"

        cd "$folder" || exit 1

        parent_drug=0
        case "$(basename "$folder" | tr '[:upper:]' '[:lower:]')" in
            abi|cbz|dox|dtx|ful|idb|lst|ptx|sn38)
                parent_drug=1
                ;;
        esac


        # 0. Completed simulation
        if [ -f MD.gro ]; then
            echo "MD.gro found. Simulation already completed."
            echo "$name DONE" >> "$ROOT/MDdone.log"
            cd "$ROOT"
            continue
        fi

        # 0b. Interrupted production MD: continue directly from checkpoint
        if [ -f MD.tpr ] && [ -f MD.cpt ]; then
            echo "MD checkpoint found. Continuing production MD."

            $GMX mdrun                 -ntmpi 1                 -ntomp "$NTOMP"                 -pin on                 -nb gpu                 -s MD.tpr                 -cpi MD.cpt                 -deffnm MD                 -v

            if [ -f MD.gro ]; then
                echo "$name DONE" >> "$ROOT/MDdone.log"
            else
                echo "$name FAILED/INCOMPLETE: MD continuation" >> "$ROOT/MDfailed.log"
            fi

            cd "$ROOT"
            continue
        fi

        # Determine the latest completed stage.
        if [ -f NPT.gro ] && [ -f NPT.cpt ]; then
            start_stage=8
            echo "NPT completed. Start from production MD."
        elif [ -f NVT.gro ] && [ -f NVT.cpt ]; then
            start_stage=7
            echo "NVT completed. Start from NPT."
        elif [ -f em.gro ]; then
            start_stage=6
            echo "EM completed. Start from NVT."
        elif [ -f sol.gro ]; then
            start_stage=5
            echo "Solvation completed. Start from index/EM."
        elif [ -f polymerdrug.gro ]; then
            start_stage=3
            echo "polymerdrug.gro found. Start from solvation."
        elif [ -f polymerdrug.pdb ]; then
            start_stage=2
            echo "polymerdrug.pdb found. Start from editconf."
        else
            start_stage=1
            echo "Starting a new simulation from Packmol."
        fi

        # 1. Packmol
        if [ "$start_stage" -le 1 ]; then
            echo "[1] Packmol"
            packmol < packmol.inp

            if [ ! -f polymerdrug.pdb ]; then
                echo "ERROR: polymerdrug.pdb was not generated."
                echo "$name FAILED: Packmol" >> "$ROOT/MDfailed.log"
                cd "$ROOT"
                continue
            fi
        fi

        # 2. PDB -> GRO
        if [ "$start_stage" -le 2 ]; then
            echo "[2] editconf"
            $GMX editconf                 -f polymerdrug.pdb                 -o polymerdrug.gro                 -box 14 14 14                 -c

            if [ ! -f polymerdrug.gro ]; then
                echo "ERROR: polymerdrug.gro was not generated."
                echo "$name FAILED: editconf" >> "$ROOT/MDfailed.log"
                cd "$ROOT"
                continue
            fi
        fi

        # 3. Solvation
        if [ "$start_stage" -le 3 ]; then
            echo "[3] solvate"
            $GMX solvate                 -cp polymerdrug.gro                 -cs spc216.gro                 -p topol.top                 -o sol.gro

            if [ ! -f sol.gro ]; then
                echo "ERROR: sol.gro was not generated."
                echo "$name FAILED: solvate" >> "$ROOT/MDfailed.log"
                cd "$ROOT"
                continue
            fi
        fi

        # 4. index.ndx
        #
        # Three grouping schemes are used:
        #   1. DTX-SI-C18+Abi-SI-C18 coassembly system
        #   2. Parent-drug systems
        #   3. Single-prodrug systems
        #
        # If a later-stage coordinate file already exists, it can be used
        # to rebuild the same index groups because atom ordering is unchanged.
        if [ ! -f index.ndx ]; then
            if [ -f sol.gro ]; then
                index_source="sol.gro"
            elif [ -f em.gro ]; then
                index_source="em.gro"
            elif [ -f NVT.gro ]; then
                index_source="NVT.gro"
            else
                index_source="NPT.gro"
            fi

            folder_name="$(basename "$folder")"

            if [ "$folder_name" = "DTX-SI-C18+Abi-SI-C18" ]; then
                echo "[4] make_ndx: DTX-SI-C18+Abi-SI-C18 coassembly system"

                printf "2 | 3 | 4 | 5 | 6 
 name 12 polymer 
 7|8 
 9|10 
 name 14 polymer 
  name 15 dtxsic18 
 name 16 othersic18 
 15|16 
 name 17 drug 
  14 | 17 
 name 18 polymer&drug 
 2|3 
 
 name 19 peg 
 
 4|5|6 
 name 20 pcl 
 
 q 
" | \
                    $GMX make_ndx -f "$index_source" -o index.ndx

            elif [ "$parent_drug" -eq 1 ]; then
                echo "[4] make_ndx: parent-drug system"

                printf "2 | 3 | 4 | 5 | 6 
 name 11 polymer 
 7|11 
  name 12 polymer&drug 
 2|3 
 
 name 13 peg 
 
 4|5|6 
 name 14 pcl 
 
 q 
" | \
                    $GMX make_ndx -f "$index_source" -o index.ndx

            else
                echo "[4] make_ndx: prodrug system"

                printf "2 | 3 | 4 | 5 | 6 
 name 12 polymer 
 7|8 
 name 13 drug 
 12 | 13 
 name 14 polymer&drug 
 2|3 
 
 name 15 peg 
 
 4|5|6 
 name 16 pcl 
 
 q 
" | \
                    $GMX make_ndx -f "$index_source" -o index.ndx
            fi
        else
            echo "[4] index.ndx already exists."
        fi

        if [ ! -f index.ndx ]; then
            echo "ERROR: index.ndx was not generated."
            echo "$name FAILED: make_ndx" >> "$ROOT/MDfailed.log"
            cd "$ROOT"
            continue
        fi


        # 5. Energy minimization
        if [ "$start_stage" -le 5 ]; then
            echo "[5] energy minimization"
            $GMX grompp                 -f em_real.mdp                 -c sol.gro                 -p topol.top                 -n index.ndx                 -o em.tpr                 -maxwarn 2 &&             $GMX mdrun -v -deffnm em

            if [ ! -f em.gro ]; then
                echo "ERROR: em.gro was not generated."
                echo "$name FAILED: EM" >> "$ROOT/MDfailed.log"
                cd "$ROOT"
                continue
            fi
        fi

        # 6. NVT
        if [ "$start_stage" -le 6 ]; then
            echo "[6] NVT"
            $GMX grompp                 -f NVT.mdp                 -c em.gro                 -p topol.top                 -n index.ndx                 -o NVT.tpr                 -maxwarn 2 &&             $GMX mdrun                 -ntmpi 1                 -ntomp "$NTOMP"                 -pin on                 -nb gpu                 -deffnm NVT                 -v

            if [ ! -f NVT.gro ] || [ ! -f NVT.cpt ]; then
                echo "ERROR: NVT did not finish correctly."
                echo "$name FAILED: NVT" >> "$ROOT/MDfailed.log"
                cd "$ROOT"
                continue
            fi
        fi

        # 7. NPT
        if [ "$start_stage" -le 7 ]; then
            echo "[7] NPT"
            $GMX grompp                 -f NPT.mdp                 -c NVT.gro                 -t NVT.cpt                 -p topol.top                 -n index.ndx                 -o NPT.tpr                 -maxwarn 2 &&             $GMX mdrun                 -ntmpi 1                 -ntomp "$NTOMP"                 -pin on                 -nb gpu                 -deffnm NPT                 -v

            if [ ! -f NPT.gro ] || [ ! -f NPT.cpt ]; then
                echo "ERROR: NPT did not finish correctly."
                echo "$name FAILED: NPT" >> "$ROOT/MDfailed.log"
                cd "$ROOT"
                continue
            fi
        fi

        # 8. Production MD
        if [ "$start_stage" -le 8 ]; then
            echo "[8] production MD"
            $GMX grompp                 -f MD.mdp                 -c NPT.gro                 -t NPT.cpt                 -p topol.top                 -n index.ndx                 -o MD.tpr                 -maxwarn 2 &&             $GMX mdrun                 -ntmpi 1                 -ntomp "$NTOMP"                 -pin on                 -nb gpu                 -deffnm MD                 -v
        fi

        if [ -f MD.gro ]; then
            echo "$name DONE" >> "$ROOT/MDdone.log"
        else
            echo "$name FAILED/INCOMPLETE: MD" >> "$ROOT/MDfailed.log"
        fi

        cd "$ROOT"

done


echo
echo "============================================================"
echo "multi_molecule_MD batch finished."
echo "Done log:   $ROOT/MDdone.log"
echo "Failed log: $ROOT/MDfailed.log"
echo "============================================================"
