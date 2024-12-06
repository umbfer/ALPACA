#! /bin/bash

scDir='/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/scripts'
scDir='/home/cattaneo/spark/PowerStatistics/scripts'

hist=hist
tempDir=tmp.$$
kmcOutputPrefix=$(mktemp ./tt.XXXXXXX)

if [ "$#" -lt 1 ]; then
    files=*.fasta
else
    files=$@
fi

minK=4
maxK=62

mkdir $tempDir
mkdir $hist

for f in $files; do
    k=$minK
    while ((k <= maxK)); do
      echo "dataset: $f k = $k"
      kmc -b -v -k$k -m2 -fm -ci0 -cs1000000 $f $kmcOutputPrefix  $tempDir
      base=$(basename $f .fasta)
      outFile=$hist/distk=${k}_${base}.hist
      kmc_dump $kmcOutputPrefix $outFile

      $scDir/hist2delta-kmc.py $outFile
      rm $outFile
      
      if (( k < 20)); then
    	  ((k+=2))
      elif (( k < 30)); then
	      ((k+=3))
      else
	      ((k+=10))
      fi
    done
done

