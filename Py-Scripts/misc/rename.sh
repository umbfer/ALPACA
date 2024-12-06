#! /bin/bash

in=/user/cattaneo/data/shigella/rnd_shigella

max=1000000

s=100

while (( s <= $max )); do
    for m in 10 25 50 75; do
	len=$(( s * m))
	SRC1=$in\_$len\_t1
	DST1=$in-T1-1000.$len.fasta
	echo "$SRC1 -- $DST1"
	hdfs dfs -mv $SRC1  $DST1

	SRC2=$in\_$len
	DST2=$in-1000.$len.fasta
	echo "$SRC2 -- $DST2"
	hdfs dfs -mv $SRC2  $DST2
    done
    s=$((s*10))
done
