#! /bin/bash

if (( $# != 3)); then
    echo "Usage: $0 seq1 seq2  hdfsDataDir"
    # exit -1
    seq1=/home/cattaneo/3rdPartiesSoftware/KMC-3.2.2/tests/ncbi_dataset/GCF_003339765.1_Mmul_1.0.fasta
    # the second sequence has been produced by makeDistance.py (same prefix -theta.fasta)
    seq2=/home/cattaneo/3rdPartiesSoftware/KMC-3.2.2/tests/ncbi_dataset/GCF_000165445.2_Mmur_3.0.fasta
    # seq2=synthetic
    dataDir=huge 
    # seq1=/home/cattaneo/3rdPartiesSoftware/KMC-3.2.2/tests/KMCCompare/S1.fasta
    # seq2=/home/cattaneo/3rdPartiesSoftware/KMC-3.2.2/tests/KMCCompare/S2.fasta
    # dataDir=data
else
    seq1=$1
    seq2=$2
    dataDir=$3
fi


cmd="spark-submit --master yarn --deploy-mode client --driver-memory 27g \
	     --num-executors 48 --executor-memory 27g --executor-cores 7 \
	     Py-Scripts/PySparkPASingleSequence.py $seq1 $seq2 $dataDir"

logFile="run-$(date '+%s').log"

echo $cmd > $logFile

$cmd >> $logFile 2>&1 


