#! /bin/bash

if (( $# != 2 )); then
    echo "Errore nei parametri. Usage:\n $0 inputDirectory outputFile"
    exit -1
fi
    
inDirectory=$1
outFile=$2

tmpFile=tt.csv

firstFile=$(ls "$inDirectory"/*.csv | head -1)

head -1 "$firstFile" > $tmpFile

for file in "$inDirectory"/*.csv; do
    grep -v "^model" "$file" >> $tmpFile
done

mv $tmpFile $outFile
