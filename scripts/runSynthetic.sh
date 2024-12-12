#! /bin/bash

export LC_NUMERIC=C

# Spark configuration variables
DRIVER_MEMORY="27g"
EXECUTOR_MEMORY="27g"
EXECUTOR_CORES=7
NUM_EXECUTORS=48

scriptDir='Py-Scripts'
# path to the directory containing input sequences
dataDir='Datasets'
# dataDir=/mnt/VolumeDati1/Dataset/PresentAbsentDatasets/ncbi_dataset

remoteDataDir=Synthetics
baseSeq=GCF_003339765.1_Mmul_1.0.fna
baseSeq='GCF_000001405.40_GRCh38.p14_coding-all.fna' # homo-sapiens coding sequence
baseSeq=fish1.fna


if (( $# < 3 )) || (( $# > 5 )); then
    echo "Usage: $0 sequence remoteDataDir local|yarn [theta [k]]"
    exit 1
fi

baseSeq=$1
remoteDataDir=$2
executionMode=$3

if [[ $executionMode != "local" && $executionMode != "yarn" ]]; then
    echo "Please specify execution mode as either 'local' or 'yarn'"
    exit -1
fi

thetaValues='0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95'
thetaValues='0.005'
kValue=""
if (( $# >= 4 )); then
    thetaValues=$4
fi
if (( $# == 5 )); then
    kValue=$5
fi



if [[ "$executionMode" == "local" ]]; then
  dataDir=$remoteDataDir
fi

seq1=${dataDir}/$baseSeq


logFile="run-$(date '+%s').log"
echo "Start Log file: $(date)e" > $logFile
echo "Log file: $logFile"

for i in $thetaValues ; do
     ${scriptDir}/makeDistance.py  ${seq1} $i
     echo "Theta: $i"
    seq2=$(printf "%s/%s-T=%.3f.fna" ${dataDir} $(basename $seq1 .fna) $i)

    if [[ "$executionMode" == "yarn" ]]; then
        cmd="spark-submit --master yarn --deploy-mode client --driver-memory $DRIVER_MEMORY \
            --num-executors $NUM_EXECUTORS --executor-memory $EXECUTOR_MEMORY --executor-cores $EXECUTOR_CORES \
            ${scriptDir}/PyPASingleSequenceOutMemory.py $seq1 $seq2 $i $remoteDataDir $executionMode $kValue"
    else
        cmd="spark-submit --master local[*] --driver-memory $DRIVER_MEMORY \
            ${scriptDir}/PyPASingleSequenceOutMemory.py $seq1 $seq2 $i $remoteDataDir $executionMode $kValue"
    fi


    echo "$(date) Comparing $seq1 vs $seq2"
    echo "$(date) Comparing $seq1 vs $seq2" >> $logFile
    echo $cmd
    $cmd >> $logFile

done


base=$(basename $baseSeq .fna)
t=0.005
tt=$(printf "%s/%s-%s-T=%.3f-T=%.3f*.csv" $dataDir $base $base $t $t)
tt=$(printf "%s/%s-%s-T=%.3f*.csv" $dataDir $base $base $t)
report=$(mktemp)

# salva l'header
head -1 $tt > $report
i=0
for f in ${dataDir}/${base}-${base}*.csv; do
    # f=$(printf "%s-%s-T=%.3f-T=%.3f*.csv" $base $base $t $t)
    echo -n "Processing file: $f -> "
    # aggiunge i risultati per ogni theta (senza header)
    tail +2 $f >> $report
    wc -l $report
    ((i++))
done

l=$(wc -l $report | cut -d ' ' -f 1)
# 8 x i + 1
tot=$((i * 8 + 1))
if ((l != tot)); then
   echo "wrong number of lines $l, expected $tot"
   wc ${dataDir}/${base}*.csv
else
    echo "$base ok $l"
fi

final="${dataDir}/${base}-$(date +%s).csv"
mv "$report" "$final"