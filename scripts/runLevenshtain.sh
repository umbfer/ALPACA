#! /bin/bash

if [ $# -ne 2 ]; then
    echo "errore numero di argomenti sbagliato $#"
    echo "Usage: runLevenshtain.sh inputFile outputDir"
    exit -1
fi

inDir=$1
outDir=$2

spark-submit --class it.unisa.di.bio.LevenshteinEditDistance \
	     --master yarn --deploy-mode client --driver-memory 16g \
	     --num-executors 48 --executor-memory 27g --executor-cores 7 \
	     --jars /home/cattaneo/.m2/repository/it/unimi/dsi/fastutil/8.5.6/fastutil-8.5.6.jar \
	     target/powerstatistics-1.0-SNAPSHOT.jar \
	     $inDir $outDir yarn


