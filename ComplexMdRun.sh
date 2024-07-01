#! /bin/bash

ARGS=1
E_BADARGS=85
E_NOFILE=86
E_Indexing=87


if [ $# -ne "$ARGS" ]
# Correct number of arguments passed to script?
then
echo "Usage: `basename $0` filename"
exit $E_BADARGS
fi
if [ ! -f "$1" ] # Does file exist?
then
echo "File \"$1\" does not exist."
exit $E_NOFILE
fi

#topology-position restraint-post-processed structure file
gmx pdb2gmx -f $1 -o conf.gro -ff gromos53a6 -water spce -ignh

#ligand gro file preparation
lline=$(cat lig.gro | wc -l) # lig.gro line 
LigAtomNum=$(expr $lline - 3) # ligand atom number
cat lig.gro | tail -n $(expr $LigAtomNum + 1) | head -n $LigAtomNum > lig.tmp
echo system | gmx genrestr -f lig.gro -o posre_lig.itp -fc 1000 1000 1000

#Receptor gro file preparation
rline=$(cat conf.gro | wc -l) # conf.gro line 
RecAtomNum=$(expr $rline - 3) # receptor atom number 
cat conf.gro | tail -n $(expr $RecAtomNum + 1) | head -n $RecAtomNum > rec.tmp

#complex gro file preparation
cat conf.gro | head -n 1 > complex.gro
ComplexNum=$(expr $LigAtomNum + $RecAtomNum)
echo " $ComplexNum" >> complex.gro
cat rec.tmp >> complex.gro
cat lig.tmp >> complex.gro
cat conf.gro | tail -n 1 >> complex.gro
rm rec.tmp lig.tmp

#Topology modification
echo "MO1                 1" >> topol.top
NumTop=$(cat topol.top | wc -l)
Num=`expr $NumTop - 21`
cat topol.top | head -n $Num > topol.tmp
echo "; Include ligand topology" >> topol.tmp
echo "#include \"lig.itp\"" >> topol.tmp
echo >> topol.tmp
echo "; Ligand position restraints" >> topol.tmp
echo "#ifdef POSRES" >> topol.tmp
echo "#include \"posre_lig.itp\"" >> topol.tmp
echo "#endif" >> topol.tmp
echo >> topol.tmp
cat topol.top | tail -n 21 >> topol.tmp
mv topol.tmp topol.top

#define the box using editconf
gmx editconf -f complex.gro -o newbox.gro -c -d 1.0 -bt dodecahedron

#Solvation
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

if [ ! -f "em.mdp" ] # Does file exist?
then
echo "File \"em.mdp\" does not exist."
exit $E_NOFILE
fi

#Adding Ions
linebond=$(grep -n bonds < topol.top | cut -d':' -f1)
linecharge=$(expr $linebond - 2)
charge=$(head -n $linecharge topol.top  | tail -n 1 | awk '{print $11}')
if [ $charge -lt  0 ]; then  atom=np; else atom=nn; fi
abscharge=$(echo $charge | tr -d -)

gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL| gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -$atom $abscharge

#Energy Minimization
if [ ! -f "em_real.mdp" ] # Does file exist?
then
echo "File \"em_real.mdp\" does not exist."
exit $E_NOFILE
fi

gmx grompp -f em_real.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em

#Indexing
echo -e "1|13\nq" | make_ndx -f em.gro -o index.ndx

ComplexIndex=$(grep Protein_MO1 index.ndx)
if [ -z "$ComplexIndex" ]
then
echo "Indexing Error, \"Protein\" and \"JZ4\" groups do not merge."
exit $E_Indexing
fi

#Equilibration-NVT
if [ ! -f "nvt.mdp" ] # Does file exist?
then
echo "File \"nvt.mdp\" does not exist."
exit $E_NOFILE
fi

gmx grompp -f nvt.mdp -c em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1
gmx mdrun -v -deffnm nvt -nt 20

#Equilibration-NPT
if [ ! -f "npt.mdp" ] # Does file exist?
then
echo "File \"npt.mdp\" does not exist."
exit $E_NOFILE
fi

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 1
gmx mdrun -v -deffnm npt -nt 20

#Production-MD
if [ ! -f "md.mdp" ] # Does file exist?
then
echo "File \"md.mdp\" does not exist."
exit $E_NOFILE
fi

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 1
gmx mdrun -v -deffnm md -nt 20
