#! /bin/bash

maxSeqNum=1000
prefix=Uniform
tmp=/tmp/aa.tmp
path=data/dataset5-1000



for len in $(seq 200000 200000 10000000); do

    for prefix in Uniform Uniform-T1 MotifRepl-U PatTransf-U; do
	    
	for gamma in 100 050 010; do

	    case $prefix in
		    Uniform | Uniform-T1)
          postfix=".fasta"
          file=$prefix-1000.$len$postfix
          one=1
          ;;
    		MotifRepl-U | PatTransf-U)
          postfix=.G=0.$gamma.fasta
          file=$prefix-1000.$len$postfix
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

		results="../k=$k/dist-k=${k}_$(basename $file .fasta).csv"
		echo "SEQ1, SEQ2, jaccard" > $results
		echo "Processing Sequences from $file for k=$k"
		for i in $(seq 1 $maxSeqNum); do
		    seqId=$(printf "$prefix-%05d" $i)
		    if [ "$one" = 1 ]; then
			seq1=$seqId-A.fasta
			seq2=$seqId-B.fasta
		    else
			seq1=$seqId.G=$gamma-A.fasta
			seq2=$seqId-G=$gamma-B.fasta
		    fi			
		    cafe -R -K $k -J /usr/local/bin/jellyfish -D jaccard -I $seq1,$seq2 > $tmp
		    line=$(tail -n 5 $tmp | head -1)
		    f1=$(echo $line | cut -d ' ' -f 1-1)
		    d=$(echo $line | cut -d ' ' -f 3-3)
		    if [ "$f1.fasta" != "$seq1" ]; then
			echo "cafe failed!!! see file $tmp"
			exit -1
		    else
			echo $seq1,$seq2, $d >> $results
		    fi
		
		done # for each pair
		echo "k = $k terminated"
	    
	    done # for each k
	    
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

