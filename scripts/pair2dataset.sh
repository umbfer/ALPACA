#!/bin/bash

if [ "$#" -eq 2 ]; then
    len=$1
    mode=$2
else
    echo "Usage:\n$0 len mode"
    exit 1
fi

baseDir=$(printf "data/dataset-%s-1000" $mode)
srcDir=$(printf "%s/len=%d" $baseDir $len)
dstDir=$baseDir/fade

hdfs dfs -mkdir $dstDir

#AltMmdels=('MotifRepl-U' 'PatTransf-U')
AltModels=('MotifRepl-Sh' 'PatTransf-Sh')
# NullModels=('Uniform' 'Uniform-T1')
NullModels=('ShuffledEColi' 'ShuffledEColi-T1')

mkdir tmp
cd tmp
for model in ${NullModels[@]}; do
    dataset=$(printf "%s/%s-????.%d.fasta" $srcDir $model $len)
    echo "get pairs from $dataset"
    hdfs dfs -get $dataset . 1>/dev/null 2>&1

    dstDataset=$(printf "%s-1000.%d.fasta" $model $len)
    localDataset=../$dstDataset
    echo "build local dataset $dstDataset"
    cat *.fasta > $localDataset

    echo "Move local dataset to $dstDir/$dstDataset"
    hdfs dfs -put $localDataset $dstDir/$dstDataset 1>/dev/null 2>&1
    rm $localDataset
    rm *.fasta
done

for model in ${AltModels[@]}; do
    for g in 010 050 100; do
	dataset=$(printf "%s/%s-????.%d.G=0.%s.fasta" $srcDir $model $len $g)
	echo "get pairs from $dataset"
	hdfs dfs -get  $dataset . 1>/dev/null 2>&1

	dstDataset=$(printf "%s-1000.%d.G=%s.fasta" $model $len $g)
	localDataset=../$dstDataset
	echo "build local dataset $dstDataset"
	cat *.fasta > $localDataset

	echo "Move local dataset to $dstDir/$dstDataset"
	hdfs dfs -put $localDataset $dstDir/$dstDataset 1>/dev/null 2>&1
	rm $localDataset
	rm *.fasta
    done
done

cd ..
rmdir tmp

hdfs dfs -ls $dstDir

echo 'Done.'
