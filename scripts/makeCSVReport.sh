#! /bin/bash

if (( $# != 3 )); then
    echo "Errore nei parametri. Usage:\n $0 inputDirectory outputFile  local|yarn"
    exit -1
fi

inDirectory=$1
outFile=$2
executionMode=$3


tmpFile=tt.csv

# shellcheck disable=SC1073
if [ $executionMode == "yarn" ]; then

  (head -1 inDirectory; grep -v "^model" inDirectory) > $tmpFile
  mv $tmpFile outFile
fi


if [ $executionMode == "local" ]; then

  firstFile=$(ls "$inDirectory"/*.csv | head -1)
  head -1 "$firstFile" > $tmpFile

  for file in "$inDirectory"/*.csv; do
      grep -v "^model" "$file" >> $tmpFile
  done

  mv $tmpFile $outFile
fi


