#! /bin/bash

if [ $# -ne 2 ]; then
    echo "errore numero di argomenti sbagliato $#"
    echo "Usage: runKMCCompare.sh kmcOutFile1 kmcOutFile2"
    exit -1
fi

inFile1=$1
inFile2=$2

spark-submit --class it.unisa.di.bio.KmerCompare \
	     --master yarn --deploy-mode client --driver-memory 16g \
	     --num-executors 48 --executor-memory 27g --executor-cores 7 \
	     target/PowerStatistics-1.3-SNAPSHOT-jar-with-dependencies.jar \
	     --mode yarn --path data/huge $inFile1 $inFile2


