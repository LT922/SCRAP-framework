#gmx mdrun -ntmpi 1 -ntomp 10 -pin on -s MD.tpr -cpi MD.cpt -deffnm MD -nb gpu -v

    gmx editconf -f polymerdrug.pdb -o polymerdrug.gro -box 14 14 14 -c

    gmx solvate -cp polymerdrug.gro -cs spc216.gro -p topol.top -o sol.gro

#    printf "2 | 3 | 4 | 5 | 6 \n name 11 polymer \n 7 \n name 12 drug \n 11 | 12 \n name 13 polymer&drug \n 2|3 \n \n name 14 peg \n \n 4|5|6 \n name 15 pcl \n \n q \n" | gmx make_ndx -f sol.gro -o index.ndx

    gmx grompp -f em_real.mdp -c sol.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 2 && gmx mdrun -v -deffnm em

    gmx grompp -f NVT.mdp -c em.gro -p topol.top -n index.ndx -o NVT.tpr -maxwarn 2 && gmx mdrun -ntmpi 1 -ntomp 10 -pin on -deffnm NVT -nb gpu -v && gmx grompp -f NPT.mdp -c NVT.gro -t NVT.cpt -p topol.top -n index.ndx -o NPT.tpr -maxwarn 2 && gmx mdrun -ntmpi 1 -ntomp 10 -pin on -deffnm NPT -nb gpu -v && gmx grompp -f MD.mdp -c NPT.gro -t NPT.cpt -p topol.top -n index.ndx -o MD.tpr -maxwarn 2 && gmx mdrun -ntmpi 1 -ntomp 10 -pin on -deffnm MD -nb gpu -v
