#! /bin/bash

tmp=tmp.$$
path=~/results/Escherichiacoli
cd $path

gammas='000 005 010 050 100 200 300 400 500'
kMin=4
kMax=33
model='EscherichiaColi'

results=$model-Jaccard.csv

declare -a mValues

k=$kMin
while ((k < kMax)); do

    for g in $gammas; do
	echo "$(date) Jaccard Distance for $model, k = $k, gamma = $g"
	jaccard.py ${model}.fasta ${model}-G=0.$g.fasta ${model}-Jaccard.csv $k

    done # for each gamma
    ((k += 2))
done # while k

rm $tmp
