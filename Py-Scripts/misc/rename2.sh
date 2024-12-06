#! /bin/bash

in=/user/cattaneo/data/shigella/

max=1000000

s=100

while (( s <= $max )); do
    for m in 10 25 50 75; do
	for g in 100 050 010; do 
	    len=$(( s * m))
	    SRC1=${in}mito-PatTransf-rn-1000.$len.G\=0.$g.fasta
	    DST1=${in}PatTransf-sh-1000.$len.G\=0.$g.fasta
	    echo "$SRC1 -- $DST1"
	    hdfs dfs -mv $SRC1  $DST1

	    SRC2=${in}mito-MotifRepl-rn-1000.$len.G\=0.$g.fasta
	    DST2=${in}MotifRepl-sh-1000.$len.G\=0.$g.fasta
	    echo "$SRC2 -- $DST2"
	    hdfs dfs -mv $SRC2  $DST2
	done
    done
    s=$((s*10))
done
