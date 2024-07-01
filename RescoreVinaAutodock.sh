for i in *out.pdbqt
do

j=`echo $i | cut -d'.' -f1`

mkdir $j
cp $i $j

cd $j

../vina_split --input $i
rm $i

for j in *.pdbqt
do 
python25 ../prepare_dpf42.py -l "$j" -r rec.pdbqt -p ga_num_evals=2500000 -p ga_pop_size=150 -p ga_run=100 -p rmstol=2.0 -p tstep=0.2 -p qstep=5.0 -p dstep=5.0 
done

cp ../rec.pdbqt .

for b in *rec.dpf 
do n=`grep -n about $b | cut -d':' -f1` 
head -n $n $b > tmp 
echo epdb >> tmp 
mv tmp $b
done

for m in ../*.map
do
cp $m .
done

cp ../rec.maps.fld .

for k in *.dpf
do
l=`echo $k | cut -d'.' -f1`
autodock4 -p $k -l $l.dlg
done

for m in *.map
do
rm $m 
done

rm rec.maps.fld

cd ..
done