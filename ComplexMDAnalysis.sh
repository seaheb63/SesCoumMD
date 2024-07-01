#! /bin/bash

#select central residue/P
gmx make_ndx -f md.tpr -n index.ndx -o index.ndx
r 
q

mkdir result 

#CONVERGENCE OF ENERGY TERMS
echo -e "Temperature\nPressure" | gmx energy -f md.edr -o Temperature-Pressure.xvg
echo -e "Total-Energy\nPotential\nKinetic-En." | gmx energy -f md.edr -o Total-Energy-Potential-Kinetic.xvg
echo Density | gmx energy -f md.edr -o Density.xvg
mv Temperature-Pressure.xvg Total-Energy-Potential-Kinetic.xvg Density.xvg result


# Remove the jumps over the boundaries
echo r_ System | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -center -pbc mol -ur compact -n index.ndx


#rerun to calculate energy
#add "energygrps = Protein MO1" to md.mdp file
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr
gmx mdrun -deffnm ie -rerun md_noPBC.xtc -nb cpu -nt 30
gmx energy -f ie.edr -o energy.xvg
???????????
mv energy.xvg result

# ANALYSIS OF RADIUS OF GYRATION
echo Protein | gmx gyrate -s md.tpr -f md_noPBC.xtc -o rg.xvg
mv rg.xvg result

# Root mean square fluctuations (RMSF)
echo Protein | gmx rmsf -f md_noPBC.xtc -s md.tpr -o rmsf.xvg -ox average.pdb -oq bfactors-residue.pdb -res -n index.ndx
mv rmsf.xvg average.pdb bfactors-residue.pdb result

# CONVERGENCE OF RMSD PROTEIN
echo Backbone Protein | gmx rms -f md_noPBC.xtc -s md.tpr -o rmsd-pro.xvg -n index.ndx
echo Backbone Backbone | gmx rms -f md_noPBC.xtc -s md.tpr -o rmsd-bb.xvg -n index.ndx
echo Backbone MO1"?-replace by number" | gmx rms -f md_noPBC.xtc -s md.tpr -o rmsd-lig.xvg -n index.ndx

mv rmsd-pro.xvg rmsd-bb.xvg rmsd-lig.xvg result


# CONVERGENCE OF RMSD COMPLEX
echo Backbone Protein_MO1 | gmx rms -f md_noPBC.xtc -s md.tpr -o rmsd-compl.xvg -n index.ndx
mv rmsd-compl.xvg result

#HYDROGEN BONDS
echo Protein Protein | gmx hbond -f md_noPBC.xtc -s md.tpr -num hb-Protein.xvg -n index.ndx
echo Protein MO1 | gmx hbond -f md_noPBC.xtc -s md.tpr -num hb-Protein-lig.xvg -n index.ndx
mv hb-Protein.xvg hb-Protein-lig.xvg result


#CLUSTER ANALYSIS
echo Backbone Protein | gmx rms -s md.tpr -f md_noPBC.xtc -m rmsd-matrix-Protein.xpm -n index.ndx -skip 50 -max 0.5
echo Backbone MO1 | gmx rms -s md.tpr -f md_noPBC.xtc -m rmsd-matrix-lig.xpm -n index.ndx -skip 50 -max 0.5
echo Backbone Protein_MO1 | gmx rms -s md.tpr -f md_noPBC.xtc -m rmsd-matrix-complex.xpm -n index.ndx -skip 50 -max 0.5

xpm2ps -f rmsd-matrix-Protein.xpm -o rmsd-matrix-Protein.eps -rainbow blue
xpm2ps -f rmsd-matrix-lig.xpm -o rmsd-matrix-lig.eps -rainbow blue
xpm2ps -f rmsd-matrix-complex.xpm -o rmsd-matrix-complex.eps -rainbow blue

echo Backbone Backbone | gmx cluster -s md.tpr -f md_noPBC.xtc -dm rmsd-matrix-Protein.xpm -dist rmsd-distribution.xvg -o clusters.xpm -sz cluster-sizes.xvg -tr cluster-transitions.xpm -ntr cluster-transitions.xvg -clid cluster-id-over-time.xvg -cl clusters.pdb -cutoff 0.2 -method gromos -n index.ndx -skip 50
echo Backbone MO1? | gmx cluster -s md.tpr -f md_noPBC.xtc -dm rmsd-matrix-lig.xpm -dist rmsd-distribution-lig.xvg -o clusters-lig.xpm -sz cluster-sizes-lig.xvg -tr cluster-transitions-lig.xpm -ntr cluster-transitions-lig.xvg -clid cluster-id-over-time-lig.xvg -cl clusters-lig.pdb -cutoff 0.2 -method gromos -n index.ndx -skip 50

mv rmsd-matrix-complex.eps rmsd-matrix-lig.eps rmsd-matrix-Protein.eps rmsd-distribution.xvg clusters.xpm  cluster-sizes.xvg cluster-transitions.xpm cluster-transitions.xvg cluster-id-over-time.xvg clusters.pdb result
mv rmsd-distribution-lig.xvg clusters-lig.xpm  cluster-sizes-lig.xvg cluster-transitions-lig.xpm cluster-transitions-lig.xvg cluster-id-over-time-lig.xvg clusters-lig.pdb result

#PCA
gmx covar -f md_noPBC.xtc -s md.tpr -ascii -n index.ndx
gmx anaeig -f md_noPBC.xtc -v eigenvec.trr -eig eigenval.xvg -s md.tpr -proj projection.xvg -first 1 -last 10 -n index.ndx
gmx anaeig -f md_noPBC.xtc -v eigenvec.trr -eig eigenval.xvg -s md.tpr -extr extreme.pdb -first 1 -last 10 -nframes 5 -n index.ndx



