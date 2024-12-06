#! /bin/bash

scriptDir='/home/cattaneo/spark/power_statistics/Py-Scripts'
dataDir='/home/cattaneo/spark/power_statistics/Datasets'
# dataDir=/mnt/VolumeDati1/Dataset/PresentAbsentDatasets/ncbi_dataset
remoteDataDir=MmulSynthetic

baseSeq=GCF_003339765.1_Mmul_1.0.fna
baseSeq=fish1.fna

seq1=${dataDir}/$baseSeq


logFile="run-$(date '+%s').log"
echo "Start Log file: $(date)e" > $logFile
echo "Log file: $logFile"

for i in 5 10 20 30 40 50 60 70 80 90 95; do

    ${scriptDir}/makeDistance.py ${seq1} $i
    seq2=${dataDir}/$(basename $seq1 .fna)-$i.fna
    
    cmd="spark-submit --master yarn --deploy-mode client --driver-memory 27g \
	     --num-executors 48 --executor-memory 27g --executor-cores 7 \
	     ${scriptDir}/PyPASingleSequenceOutMemory.py $seq1 synthetic $i $remoteDataDir"

    
    echo "Comparing $(basename $seq1) vs $(basename $seq2)"
    echo "Comparing $seq1 vs $seq2" >> $logFile
    echo $cmd >> $logFile

    $cmd >> $logFile 2>&1 
done

