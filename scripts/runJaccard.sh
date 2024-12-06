#! /bin/bash

maxSeqNum=1000
cores=8
prefix=""
tmp=/tmp/aa.$$
path=data/dataset5-1000
resultsPath='/home/cattaneo/spark/power_statistics/data/results/Jaccard'


for len in $(seq 4800000 200000 10000000); do

    for prefix in Uniform Uniform-T1 MotifRepl PatTransf; do
	    
	for gamma in 100 050 010; do

	    case $prefix in
		Uniform | Uniform-T1)
		    postfix=".fasta"
		    file=$prefix-1000.$len$postfix
		    one=1
		    ;;
		MotifRepl | PatTransf)
		    postfix=.G=0.$gamma.fasta
		    file=$prefix-U-1000.$len$postfix
		    one=0
		    ;;
		*)
		    echo "Unknown prefix"
		    exit -1
		    ;;
	    esac

	    echo "Downloading file: <$path/$file>"
	    hdfs dfs -get $path/$file . > /dev/null 2>&1
	    
	    echo "Splitting file"
	    splitFasta.py $file > /dev/null
		
	    cd seqDists
	    
	    for k in 4 6 8 10; do

		results=${resultsPath}"/k=$k/dist-k=${k}_$(basename $file .fasta).csv"
		echo "SEQ1, SEQ2, jaccard" > $results
		results2=$(dirname $results)/$(basename $results .csv)-values.csv
		echo "A, B, C, D, N, N-D-A" > $results2
		echo "Processing Sequences from $file for k=$k"
		
		for i in $(seq 1 $maxSeqNum); do
		    seqId=$(printf "$prefix-%05d" $i)
		    if [ "$one" = 1 ]; then
			seq1=$seqId-A.fasta
			seq2=$seqId-B.fasta
		    else
			seq1=$seqId.G=$gamma-A.fasta
			seq2=$seqId.G=$gamma-B.fasta
		    fi

		    runCafe.sh $seq1 $seq2 $k $results & # run in background

		    bck=( $(jobs -p) )
		    while (( ${#bck[@]} >= $cores )); do
			wait -n
			bck=( $(jobs -p) )
		    done
		
		done # for each pair
		echo "k = $k terminated"
	    
	    done # for each k
	    # wait for all background processes
	    bck=( $(jobs -p) )
	    while (( ${#bck[@]} > 0 )); do
		wait -n
		bck=( $(jobs -p) )
	    done
	    echo "$file processed."
	    cd ..
	    rm -fr seqDists
	    rm $file
	    
	    if [ "$one" = "1" ]; then
		break;
	    fi
	
	done # for each gamma
    
    done # for each model

done #for each len


