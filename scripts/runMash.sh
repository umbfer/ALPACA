#! /bin/bash

tmp=tmp.$$
path=~/results/Escherichiacoli
cd $path

gammas='000 005 010 050 100 200 300 400 500'
sketches='1000 10000 100000'
kMin=4
kMax=33
model='EscherichiaColi'

results=$model-Mash.csv
if [ ! -s $results ]; then
    echo "model,gamma,k,sketchSize,MashDistance,p-value, A/sketchsize" > $results
fi

declare -a mValues

for sketchSize in $sketches; do
    k=$kMin
    while ((k < kMax)); do
	echo "$(date) Sketch for k = $k size = $sketchSize"
	mash sketch -s $sketchSize -k $k ${model}.fasta
	
	for g in $gammas; do
	    echo "Distances for gamma = $g"
	    mash sketch -s $sketchSize -k $k ${model}-G=0.$g.fasta
	    mash dist ${model}.fasta.msh ${model}-G\=0.$g.fasta.msh > $tmp
	    mValues=($(cat $tmp | cut -f 3-5))
	    distance=${mValues[0]}
	    pvalue=${mValues[1]}
	    A=${mValues[2]}
	    echo "$model, 0.$g, $k, $sketchSize, $distance, $pvalue, $A" >> $results
	    
	done # for each gamma
	((k += 2))
    done # while k
done # for each sketch
rm $tmp
