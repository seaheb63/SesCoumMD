#! /bin/bash

#ligand & receptor preparation

ARGS=2
E_BADARGS=85
E_NOFILE=86


if [ $# -ne "$ARGS" ] # Correct number of arguments passed to script?
then
echo "Usage: `basename $0` receptor ligand"
exit $E_BADARGS
fi
if [ ! -f "$1" ] # Does file exist?
then
echo "File \"$1\" does not exist."
exit $E_NOFILE
fi
echo
echo "***************************************************************************************************"
echo "***     Crystalloghraphic ligand in pdb file ($2) will be used to determine GRID BOX center     ***"
echo "***************************************************************************************************"

#CORE NUMBER
echo "Enter CPU Core Number to Use:"
read cnum

#receptor preparation
prepare_receptor4.py -r $1 -o $1"qt" -A checkhydrogens

#ligand preparation
Dir=$(pwd)
cd lig #ligand pdb/xyz file sould be in lig directory

for i in *; do
j=$(echo $i | cut -d'.' -f2)
k=$(echo $i | cut -d'.' -f1)
if [ $j == "xyz" ]
then
   obabel "$k".xyz -O "$k".mol2 -p7.4
elif [ $j == "pdb" ]
then 
   obabel "$k".pdb -O "$k".mol2 -p7.4
else
   echo "********NO LIGAND********"
fi

#LIG PDBQT
prepare_ligand4.py -l "$k".mol2

##DPF
prepare_dpf42.py -l "$k".pdbqt -r $1"qt" -o $k.dpf -p ga_num_evals=2500000 -p ga_pop_size=150 -p ga_run=100 -p rmstol=2.0 -p tstep=0.2 -p qstep=5.0 -p dstep=5.0

mv $k.pdbqt $Dir
mv $k.dpf $Dir

done 
cd $Dir

#GPF
prepare_gpf4.py -l $2"qt" -r $1"qt" -p ligand_types='A,Br,Cl,C,F,HD,H,HS,I' -y -o g1.gpf 
echo "npts 60 60 60                        # num.grid points in xyz" > g1.tmp
g1line=$(grep -n dielectric < g1.gpf | cut -d':' -f1)
tail -n $(expr $g1line - 1) g1.gpf >> g1.tmp
mv g1.tmp g1.gpf

prepare_gpf4.py -l $2"qt" -r $1"qt" -p ligand_types='NA,N,NS,OA,OS,P,SA,S' -y -o g2.gpf
echo "npts 60 60 60                        # num.grid points in xyz" > g2.tmp
g2line=$(grep -n dielectric < g2.gpf | cut -d':' -f1)
tail -n $(expr $g2line - 1) g2.gpf >> g2.tmp
mv g2.tmp g2.gpf

#GRID MAP
autogrid4 -p g1.gpf -l g.glg
autogrid4 -p g2.gpf -l g.glg

#PARALLEL RUN
TaskLimit=$cnum
Task=0
for i in *.dpf; do
{
j=`echo $i | cut -d'.' -f1` 
autodock4 -p $i -l $j.dlg &
Task=$(($Task+1))
if [ $Task -ge $TaskLimit ]
then 
wait
Task=0
fi 
}
done


