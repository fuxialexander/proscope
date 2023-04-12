# a shell script to run a molecular dynamics simulation using gromacs
# the structure files are in the current directory and is a multi-chain pdb file
# generate the topology file using pdb2gmx
# the topology file is named topol.top
# the structure file is named structure.pdb
# the force field is a99SB-disp
structure_pdb=$1
structure=$(basename $structure_pdb .pdb)
source /usr/local/gromacs/bin/GMXRC
# run pdb2gmx, ensure the gro file keeps the chain information
rm \#*
gmx pdb2gmx -f $structure.pdb -ignh -o $structure.ff.pdb -p $structure.topol.top -ff a99SBdisp -water select << INPUT
1
INPUT

# define the box
echo "Define the box: Box size is 1 nm"
gmx editconf -f $structure.ff.pdb -o $structure.ff_box.pdb -bt cubic -d 1.0 

# solvate the box
echo "Solvate the box: Solvent is TIP4P" 
gmx solvate -cp $structure.ff_box.pdb -cs tip4p.gro -o $structure.ff_box_solv.pdb -p $structure.topol.top

# generate the ions
echo "Generate the ions: Na+ and Cl-"
gmx grompp -f ions.mdp -c $structure.ff_box_solv.pdb -p $structure.topol.top -o $structure.ions.tpr -maxwarn 1
gmx genion -s $structure.ions.tpr -o $structure.ff_box_solv_ions.pdb -p $structure.topol.top -pname NA -nname CL -neutral << INPUT
13
INPUT

# Assemble the binary input using grompp using minim.mdp
echo "Minimization"
gmx grompp -f em.mdp -c $structure.ff_box_solv_ions.pdb -p $structure.topol.top -o $structure.em.tpr

# run the minimization
gmx mdrun -v -deffnm $structure.em -nt 16 -nb gpu -c $structure.em.pdb -e $structure.em.edr -g $structure.em.log
gmx energy -f $structure.em.edr -o $structure.potential.xvg << INPUT
10 0
INPUT


# NVT
echo "NVT"
gmx grompp -f nvt.mdp -c $structure.em.pdb -r $structure.em.pdb -p $structure.topol.top -o $structure.nvt.tpr
gmx mdrun -v -deffnm $structure.nvt -nt 16 -nb gpu -c $structure.nvt.pdb -e $structure.nvt.edr -g $structure.nvt.log
gmx energy -f $structure.nvt.edr -o $structure.temperature.xvg << INPUT
16 0
INPUT

# NPT
echo "NPT"
gmx grompp -f npt.mdp -c $structure.nvt.pdb -r $structure.nvt.pdb -t $structure.nvt.cpt -p $structure.topol.top -o $structure.npt.tpr
gmx mdrun -v -deffnm $structure.npt -nt 16 -nb gpu -c $structure.npt.pdb -e $structure.npt.edr -g $structure.npt.log
gmx energy -f $structure.npt.edr -o $structure.pressure.xvg << INPUT
18 0
INPUT
gmx energy -f $structure.npt.edr -o $structure.density.xvg << INPUT
24 0
INPUT

# MD
echo "50ns MD"
gmx grompp -f md.mdp -c $structure.npt.pdb -t $structure.npt.cpt -p $structure.topol.top -o $structure.md.tpr
gmx mdrun -v -deffnm $structure.md -nt 16 -nb gpu -c $structure.md.pdb -e $structure.md.edr -g $structure.md.log -x $structure.md.xtc -o $structure.trr

