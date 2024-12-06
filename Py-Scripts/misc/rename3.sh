#! /bin/bash


k=11

for len in 5000000 7500000 10000000 ; do
    for am in MotifRepl PatTransf; do
	for g in 010 050 100; do
	    SRC1=dist-k\=${k}_mito-${am}-rn-1000.${len}.G=0.${g}.csv
	    DST1=dist-k\=${k}_${am}-mi-1000.${len}.G=0.${g}.csv
	    echo "$SRC1 -- $DST1"
	    mv $SRC1  $DST1
	done
    done
done
