#! /bin/bash
eval "$(conda shell.bash hook)"
conda activate drugfilter_py39_env

omp=30
num=`sed -n '$=' $1`
avgnum=$[(num+omp-1)/$omp]
filenum=$[(num+avgnum-1)/$avgnum]

split --numeric=1 -l $avgnum $1 liglist_
ls|grep liglist_|xargs -n1 -i{} mv {} {}.smi

if [ $filenum -le 1 ]; then
#cat liglist_01.smi | tail -n +1 >> Cal_ademt.smi
python ./src/proxy.drug.admet.fix6.py liglist_01.smi Cal_ademt.csv
else
#filenum=`printf "%02d\n" $filenum`
#for i in `eval echo {01..$filenum}`;do echo $i; done
#parallel -j $filenum "cat liglist_{}.smi | tail -n +1 >> Cal_ademt.smi" ::: `eval echo {01..$filenum}`
parallel -j $filenum "python ./src/proxy.drug.admet.fix8.py liglist_{}.smi liglist_{}.csv" ::: `eval echo {01..$filenum}`
cat liglist_01.csv | head -n 1 >> Cal_ademt.csv
parallel -j $filenum "cat liglist_{}.csv | tail -n +2 >> Cal_ademt.csv" ::: `eval echo {01..$filenum}`
fi

rm -rf liglist*