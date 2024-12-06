#! /usr/bin/bash


base=$1
t=5
tt=$(printf "%s-%s-%02d-%02d*.csv" $base $base $t $t)
report=${base}-$(date +%s).csv

# salva l'header
head -1 $tt > $report
i=0
for f in ${base}-${base}*.csv; do
    # f=$(printf "%s-%s-%02d-%02d*.csv" $base $base $t $t)
    echo -n "Processing file: $f -> "
    # aggiunge i risoltati per ogni theta (senza header)
    tail +2 $f >> $report
    wc -l $report
    ((i++))
done

l=$(wc -l $report | cut -d ' ' -f 1)
# 8 x i + 1
tot=$((i * 8 + 1))
if (($l != $tot)); then
   echo "wrong number of lines $l"
   wc ${base}*.csv
else
    echo $base ok $l
fi
