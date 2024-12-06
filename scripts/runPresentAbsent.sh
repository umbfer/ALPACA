#! /bin/bash

# Spark configuration variables
DRIVER_MEMORY="27g"
EXECUTOR_MEMORY="27g"
EXECUTOR_CORES=7
NUM_EXECUTORS=48

# Check input arguments
if (( $# != 3 )); then
    echo "Usage: $0 seqLen mainDataDir local|yarn"
    exit -1
else
    seqLen=$1
    dataDir=$2
    executionMode=$3

    if [[ $executionMode != "local" && $executionMode != "yarn" ]]; then
        echo "Please specify execution mode as either 'local' or 'yarn'"
        exit -1
    fi
fi

# Base Spark command
cmd="spark-submit --master $executionMode --deploy-mode client --driver-memory $DRIVER_MEMORY"

# Add executor configuration for yarn mode
if [[ $executionMode == "yarn" ]]; then
    cmd="$cmd --num-executors $NUM_EXECUTORS --executor-memory $EXECUTOR_MEMORY --executor-cores $EXECUTOR_CORES"
fi


# Add script and parameters
cmd="$cmd Py-Scripts/PySparkPresentAbsent4.py $seqLen $dataDir"

echo $cmd
# Logging
logFile="run-${seqLen}-$(date '+%s').log"
echo $cmd > $logFile

# Execute command
$cmd >> $logFile 2>&1
