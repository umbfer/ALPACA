#! /bin/bash

inFile=$1
len=2000000

if test ! -f $inFile; then
    echo "Cannot find input file: $inFile"
    exit -1
fi
winSize=100
i=1
j=1
while ((i < ((len - i)) && ((j <= 10)))) ; do
    echo $i
    # salta la prima riga con l'intestazione
    awk ' NR==2 { print $1}' $inFile | tail -c +$i   | head -c $winSize > $(basename $inFile .fasta)-${j}.fasta
    ((i+=winSize))
    ((j++))
done
