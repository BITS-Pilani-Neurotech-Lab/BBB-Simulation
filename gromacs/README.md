# GROMACS Commands Reference

No `.mdp` parameter files were found in the original project folder, so
`gromacs/mdp/` is currently empty — add the ones actually used for this
project's minimization/equilibration/production runs. Typical GROMACS
workflow commands for reference once they're in place:

```bash
# Ion placement (neutralize + add salt)
gmx grompp -f mdp/ions.mdp -c system.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o system_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Energy minimization
gmx grompp -f mdp/minim.mdp -c system_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

# NVT equilibration
gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT equilibration
gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Production MD
gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md

# Post-processing: remove PBC jumps, center the protein
gmx trjconv -s md.tpr -f md.xtc -o md_centered.xtc -pbc mol -center
```

Adjust flags (e.g. `-ntmpi`, `-ntomp`, `-gpu_id`) to match your hardware.
