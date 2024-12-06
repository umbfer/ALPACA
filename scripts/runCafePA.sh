#! /bin/bash

scDir='/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/scripts'
scDir='/home/cattaneo/spark/PowerStatistics/scripts'

hist=hist
#tempDir=tmp.$$
jfOutput=$(mktemp ./tt.XXXXXXX)
res1=distances.txt
res2=jaccardValues.txt

if [ "$#" -lt 1 ]; then
    models='Uniform MotifRepl PatTransf Uniform-T1'
else
    models=$@
fi

lens='2000 20000 200000 2000000 20000000'

gVals='10 50 100'

minK=4
maxK=62


for m in $models; do
  case $m in
	  "Uniform*")
	    gammas=''
	    ;;
	  "MotifRepl" | "PatTransf")
	    gammas=$gVals
	    ;;
  esac
	    
  for len in $lens; do
    for seqId in $(seq 1 20); do
      for g in $gammas; do
        if [ -z "$g" ]; then
            gamma=''
        else
            gamma=$(printf ".G=%03d" $g)
        fi
        k=$minK
        while ((k <= maxK)); do
            dataset=$(printf %s-%04d.%d%s $m $seqId $len $gamma)
            seq1=$dataset-A.fasta
            seq2=$dataset-B.fasta

            cmd="cafe-mod -R -K $k -J /usr/local/bin/jellyfish -D jaccard -I $seq1,$seq2"
            echo "Processing $cmd"

            $cmd > $jfOutput

            echo "Command finished"
            rm hash_${dataset}*
            line=$(tail -n 5 $jfOutput| head -1)
            f1=$(echo $line | cut -d ' ' -f 1-1)
            d=$(echo $line | cut -d ' ' -f 3-3)
            if [ "$f1.fasta" != $(basename "$seq1") ]; then
              echo "cafe failed!!! see file $jfOutput"
              exit -1
                else
              echo $seq1,$seq2, $d >> $res1
              line2=$(tail -n 9 $jfOutput | head -1)
              echo $line2 >> $res2
            fi

            if (( k < 20)); then
              ((k+=2))
            elif (( k < 30)); then
              ((k+=3))
            else
              ((k+=10))
            fi
        done
      done
    done
  done
done

