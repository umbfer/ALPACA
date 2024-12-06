#! /usr/bin/python3

import itertools
import subprocess
import time


scriptPath = '/home/cattaneo/spark/power_statistics/Py-Scripts/PyPASingleSequenceOutMemory.py'
# dataDir='/home/cattaneo/spark/power_statistics/Dataset'
dataDir = '/mnt/VolumeDati1/Dataset/PresentAbsentDatasets/ncbi_dataset'
remoteDataDir = 'huge'


seqs = [ 'GCF_000165445.2_Mmur_3.0.fasta', 
	 'GCF_000955945.1_Caty_1.0.fasta',
	 'GCF_003339765.1_Mmul_1.0.fasta',
	 'GCF_012559485.2_MFA.2.fasta']

seqs = [ 'fish1.fasta', 'fish2.fasta'] 


for p in itertools.combinations( seqs, 2):
	 
    seq1 = f"{dataDir}/{p[0]}"
    seq2 = f"{dataDir}/{p[1]}"

    print(f"Running {p[0]} vs {p[1]}")
    
    cmd = f"spark-submit --master yarn --deploy-mode client --driver-memory 27g \
	     --num-executors 48 --executor-memory 27g --executor-cores 7 \
	     {scriptPath} {seq1} {seq2} {remoteDataDir}"

    logFile = f"run-{int(time.time())}.log"
    # $cmd >> $logFile 2>&1 
    out_f = open(logFile, 'w')
    subprocess.run( cmd.split(), stdout = out_f, text = True, stderr = subprocess.STDOUT)



