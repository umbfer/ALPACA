#! /bin/bash

in=/user/cattaneo/data/shigella/rnd_shigella

max=1000000

s=100

for f in $(ls *mito*.csv); do
    n=$(echo $f | sed 's/rnd_mito/mitocondri_q10/')

    echo "$f -> $n"
    mv $f $n
done
