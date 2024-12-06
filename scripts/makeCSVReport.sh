#! /bin/bash

if (( $# != 2 )); then
    echo "Errore nei parametri. Usage:\n $0 inputDirectory outputFile"
    exit -1
fi

inDirectory=$1
outFile=$2

tmpFile=tt.csv

(head -1 inDirectory; grep -v "^model" inDirectory) > $tmpFile

mv $tmpFile outFile
