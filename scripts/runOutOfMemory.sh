#! /bin/bash

dataDir='/home/cattaneo/spark/power_statistics/Dataset'
dataDir=/mnt/VolumeDati1/Dataset/PresentAbsentDatasets/ncbi_dataset
dataDir=/mnt/VolumeDati1/Dataset/PresentAbsentDatasets/tests
# GCF_000165445.2_Mmur_3.0.fasta
# GCF_000955945.1_Caty_1.0.fasta
# GCF_003339765.1_Mmul_1.0.fasta
# GCF_012559485.2_MFA.2.fasta


if (( $# != 3)); then
    seq1=${dataDir}/GCF_000165445.2_Mmur_3.0.fasta
    seq1=${dataDir}/S1.fasta
    # synthetic => the second sequence has been produced by makeDistance.py (same prefix -theta.fasta)
    # seq2=synthetic
    seq2=${dataDir}/GCF_000955945.1_Caty_1.0.fasta
    seq2=${dataDir}/S2.fasta
    remoteDataDir=huge 
    # dataDir=tests
    echo "Usage: $0 seq1 seq2  hdfsDataDir"
    echo "Using: $0 $seq1 $seq2 $remoteDataDir"
    # exit -1
else
    seq1=${dataDir}/$1
    seq2=${dataDir}/$2
    remoteDataDir=$3
fi


cmd="spark-submit --master yarn --deploy-mode client --driver-memory 27g \
	     --num-executors 48 --executor-memory 27g --executor-cores 7 \
	     Py-Scripts/PyPASingleSequenceOutMemory.py $seq1 $seq2 $remoteDataDir"

logFile="run-$(date '+%s').log"

echo $cmd > $logFile

$cmd >> $logFile 2>&1 


