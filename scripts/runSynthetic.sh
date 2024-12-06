#! /bin/bash

scriptDir='/home/cattaneo/spark/power_statistics/Py-Scripts'
dataDir='/home/cattaneo/spark/power_statistics/Datasets'
# dataDir=/mnt/VolumeDati1/Dataset/PresentAbsentDatasets/ncbi_dataset
remoteDataDir=Synthetics
baseSeq=GCF_003339765.1_Mmul_1.0.fna
baseSeq='GCF_000001405.40_GRCh38.p14_coding-all.fna' # homo-sapiens coding sequence
baseSeq=fish1.fna


if (( $# < 2)) || (($# > 4)); then
    echo "Usage: $0 sequence remoteDataDir [theta [k]]"
    exit -1
else
    if (($# >= 2)) ; then
	baseSeq=$1
	remoteDataDir=$2
	thetaValues='5 10 20 30 40 50 60 70 80 90 95'
	kValue=""
    fi
    if (($# >= 3)) ; then
	thetaValues=$3
    fi
    if (($# >= 4)) ; then
	kValue=$4
    fi
fi

seq1=${dataDir}/$baseSeq


logFile="run-$(date '+%s').log"
echo "Start Log file: $(date)e" > $logFile
echo "Log file: $logFile"

for i in $thetaValues ; do

    /usr/local/bin/MoveAway  ${seq1} $i
    seq2=$(printf "%s/%s-%02d.fna" ${dataDir} $(basename $seq1 .fna) $i)
    
    cmd="spark-submit --master yarn --deploy-mode client --driver-memory 27g \
	     --num-executors 48 --executor-memory 27g --executor-cores 7 \
	     ${scriptDir}/PyPASingleSequenceOutMemory.py $seq1 $seq2 $i $remoteDataDir $kValue"
    
    echo "$(date) Comparing $seq1 vs $seq2"
    echo "$(date) Comparing $seq1 vs $seq2" >> $logFile
    $cmd >> $logFile

done

