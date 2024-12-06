#! /usr/local/bin/python3

import re
import os
import os.path
from pathlib import Path
import sys
import glob
import tempfile
import shutil
import copy
import subprocess
import math
import csv
import time
import makeDistance as mkd

import numpy as np
import py_kmc_api as kmc

from operator import add
import pyspark
from pyspark.sql import SparkSession
from pyspark import SparkFiles
from pyspark.sql.types import StructType, StructField, StringType, IntegerType, LongType
import pyspark.sql.functions as sf


hdfsPrefixPath = 'hdfs://master2:9000/user/cattaneo'
# hdfsPrefixPath = '/Users/pipp8'

hdfsDataDir = ''
spark = []
sc = []
thetaValue = 0

nTests = 1000
minK = 4
maxK = 32
stepK = 4
sketchSizes = [1000, 10000, 100000]
# sketchSizes = [10000]

outFilePrefix = 'PresentAbsentRealGenomeData'

# global broadcast variables
totDistinctKmerAAcc = []
totDistinctKmerBAcc = []
totKmerAAcc = []
totKmerBAcc = []
kmerStats = []



class EntropyData:
    def __init__(self, nKeys, totalKmerCnt, Hk):
        self.nKeys = nKeys
        self.totalKmerCnt = totalKmerCnt;
        self.Hk = Hk

    def getDelta(self):
        return float(self.nKeys) / (2 * self.totalKmerCnt)

    def getError(self):
        try:
            result = self.getDelta() / self.Hk
        except ZeroDivisionError:
            result = 0
        return result

    def toString(self):
        return [self.nKeys, 2 * self.totalKmerCnt, self.getDelta(), self.Hk, self.getError()]




class MashData:
    def __init__(self, cmdResults):
        mr = cmdResults.split()
        mashAN = mr[4].decode('UTF-8')
        self.Pv = float(mr[2])
        self.dist = float(mr[3])
        try:
            ns = mashAN.index('/')
            self.A = int(mashAN[:ns])
            self.N = int(mashAN[ns+1:])
        except ValueError:
            self.A = 0
            self.N = 0

    def toString( self):
        return [ self.Pv, self.dist, self.A, self.N ]


def checkPathExists(path: str) -> bool:
    global hdfsDataDir, spark
    # spark is a SparkSession
    sc = spark.sparkContext
    fs = sc._jvm.org.apache.hadoop.fs.FileSystem.get(
        sc._jvm.java.net.URI.create(hdfsDataDir),
        sc._jsc.hadoopConfiguration(),)
    return fs.exists(sc._jvm.org.apache.hadoop.fs.Path(path))


def hamming_distance(seq1: str, seq2: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def hamming_distance2(seq1: str, seq2: str) -> int:
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(seq1, seq2))))




# calcola i valori dell'entropia per non caricare due volte l'istogramma




# we use numpy to not reimplment z-score stndardization from scratch
def NormalizedSquaredEuclideanDistance( vector):
    # (tot1, tot2) = (0, 0)
    # for x in vector:
    #     tot1 += x[0]
    #     tot2 += x[1]
    # n = len(vector)
    # mean1 = tot1 / n
    # mean2 = tot2 / n
    # (totDifferences1,totDifferences2) = (0,0)
    # for v in [((value[0] - mean1)**2, (value[1] - mean2)**2)  for value in vector]:
    #     totDifferences1 += v[0]
    #     totDifferences2 += v[1]
    # standardDeviation1 = (totDifferences1 / n) ** 0.5
    # standardDeviation2 = (totDifferences2 / n) ** 0.5
    # zscores = [((v[0] - mean1) / standardDeviation1, (v[1] - mean2) / standardDeviation2) for v in vector]

    # avg = np.mean( vector, axis=1)
    # std = np.std( vector, axis=1)
    #
    # z0_np = (vector[0] - avg[0]) / std[0]
    # z1_np = (vector[1] - avg[1]) / std[1]
    # tot = 0
    # for i in range(vector.shape[1]):
    #     tot += ((z0_np[i] - z1_np[i]) ** 2)
    #
    # ZEu = tot ** 0.5
    #
    # m = vector.shape[1]
    # D = 2 * m * (1 - (np.dot(vector[0], vector[1]) - m * avg[0] * avg[1]) / (m * std[0] * std[1]))
    var = np.var( vector, axis=1)

    NED = 0.5 * np.var(vector[0] - vector[1]) / (var[0] + var[1])
    return NED




# run jaccard on sequence pair ds with kmer of length = k
def runPresentAbsent(  bothCnt: int, leftCnt: int, rightCnt: int, k: int):

    A = int(bothCnt)
    B = int(leftCnt)
    C = int(rightCnt)

    NMax = pow(4, k)
    M01M10 = leftCnt + rightCnt
    M01M10M11 = bothCnt + M01M10
    absentCnt = NMax - (A + B + C) # NMax - M01M10M11
    D = absentCnt
    # (M10 + M01) / (M11 + M10 + M01)

    # Anderberg dissimilarity => Anderberg = 1 - (A/(A + B) + A/(A + C) + D/(C + D) + D/(B + D))/4
    try:
        anderberg = 1 - (A/float(A + B) + A/float(A + C) + D/float(C + D) + D/float(B + D))/4.0
    except (ZeroDivisionError, ValueError):
        anderberg = 1.000001

    # Antidice dissimilarity => Antidice = 1 - A/(A + 2(B + C))
    try:
        antidice = 1 - A / float(A + 2.0 * (B + C))
    except (ZeroDivisionError, ValueError):
        antidice = 1.000001

    # Dice dissimilarity => Dice = 1 - 2A/(2A + B + C)
    try:
        dice = 1 - 2*A / float(2.0*A + B + C)
    except (ZeroDivisionError, ValueError):
        dice  = 1.000001
    # Gower dissimilarity => Gower = 1 - A x D/sqrt(A + B) x(A + C) x (D + B x (D + C)
    try:
        gower = 1 - A * D / math.sqrt((A + B) * (A + C) * (D + B * (D + C)))
    except (ZeroDivisionError, ValueError):
        gower = 1.000001

    # Hamman dissimilarity => Hamman = 1 - [((A + D) - (B + C))/N]2
    try:
        hamman = 1 - math.pow((((A + D) - (B + C)) / float(NMax)), 2.0)
    except (ZeroDivisionError, ValueError):
        hamman = 1.000001

    # Hamming dissimilarity => Hamming = (B + C)/N
    try:
        hamming = (B + C)/ NMax
    except (ZeroDivisionError, ValueError):
        hamming = 1.000001

    # Jaccard dissimilarity => Jaccard = 1 - A/(N - D)
    try:
        jaccard = 1 - A / (NMax - D)
    except (ZeroDivisionError, ValueError):
        jaccard = 1.000001

    try:
        jaccardDistance = 1 - min( 1.0, A / float(NMax - D))
    except (ZeroDivisionError, ValueError):
        jaccardDistance = 1.000001


    # Kulczynski dissimilarity => Kulczynski = 1 - (A/(A + B) + A/(A + C)) / 2
    try:
        kulczynski = 1 - (A / float(A + B) + A / float(A + C)) / 2.0
    except (ZeroDivisionError, ValueError):
        kulczynski = 1.000001

    # Matching dissimilarity => Matching = 1 - (A + D)/N
    try:
        matching = 1 - (A + D) / NMax
    except (ZeroDivisionError,ValueError):
        matching = 1.000001

    # Ochiai dissimilarity => Ochiai = 1 - A/sqrt(A + B) x (A + C)
    try:
        ochiai = 1 - A / math.sqrt((A + B) * (A + C))
    except (ZeroDivisionError, ValueError):
        ochiai = 1.000001

    # Phi dissimilarity => Phi = 1 - [(A x  B x  C x D)/sqrt(A + B) x (A + C) x (D + B) x (D + C)]2
    try:
        phi = 1 - math.pow((A * B * C * D)/ math.sqrt((A + B) * (A + C) * (D + B) * (D + C)), 2.0)
    except (ZeroDivisionError, ValueError):
        phi = 1.000001

    # Russel dissimilarity => Russel = 1 - A/N
    try:
        russel = 1 - A / NMax
    except (ZeroDivisionError, ValueError):
        russel = 1.000001

    # Sneath dissimilarity => Sneath = 1 - 2(A + D)/(2 x (A + D) + (B + C))
    try:
        sneath = 1 - 2.0 * (A + D) / (2.0 * (A + D) + (B + C))
    except (ZeroDivisionError, ValueError):
        sneath = 1.000001

    # Tanimoto dissimilarity => Tanimoto = 1 - (A + D)/((A + D) + 2(B + C))
    try:
        tanimoto = 1 - (A + D) / float((A + D) + 2.0 * (B + C))
    except (ZeroDivisionError, ValueError):
        tanimoto = 1.000001

    # Yule dissimilarity => Yule = 1 - [(A x D - B x C)/(A x D + B x C)]2
    try:
        yule = 1 - math.pow(((A * D - B * C) / float(A * D + B * C)), 2.0)
    except (ZeroDivisionError, ValueError):
        yule = 1.000001

    # salva il risultato nel file CSV
    # dati present / absent e distanze present absent
    data1 = [ A, B, C, str(D), str(NMax), str(A/NMax),
             anderberg, antidice, dice, gower, hamman, hamming, jaccard,
             kulczynski, matching, ochiai, phi, russel, sneath, tanimoto, yule]

    return data1




def runMash(inputDS1, inputDS2, k):
    # run mash on the same sequence pair
    mashValues = []
    for i in range(len(sketchSizes)):
        # extract mash sketch from the first sequence
        cmd = "/usr/local/bin/mash sketch -s %d -k %d %s" % (sketchSizes[i], k, inputDS1)
        p = subprocess.Popen(cmd.split())
        p.wait()

        cmd = "/usr/local/bin/mash sketch -s %d -k %d %s" % (sketchSizes[i], k, inputDS2)
        p = subprocess.Popen(cmd.split())
        p.wait()

        cmd = "/usr/local/bin/mash dist %s.msh %s.msh" % (inputDS1, inputDS2)
        out = subprocess.check_output(cmd.split())

        mashValues.append( MashData( out))

    # dati mash distance
    data2 = []
    for i in range(len(sketchSizes)):
        data2 = data2 + mashValues[i].toString()

    # clean up remove kmc temporary files
    os.remove(inputDS1 + '.msh')
    os.remove(inputDS2 + '.msh')

    return data2




# load histogram for both sequences (for counter based measures such as D2)
# and calculate Entropy of the sequence
# dest file è la path sull'HDFS già nel formato hdfs://host:port/xxx/yyy
def loadHistogramOnHDFS(histFile: str, destFile: str):

    # tmp = histFile + '.txt'
    #
    # # dump the result -> kmer histogram
    # cmd = "/usr/local/bin/kmc_dump %s %s" % (histFile, tmp)
    # p = subprocess.Popen(cmd.split())
    # p.wait()
    # print(f"****** Dumping {histFile} kmers counting ******")
    #
    # # trasferisce sull'HDFS il file testuale
    # print(f"****** Transferring to hdfs {histFile} -> {destFile} ******")
    # cmd = "hdfs dfs -put %s %s" % (tmp, destFile)
    # p = subprocess.Popen(cmd.split())
    # p.wait()
    # os.remove(tmp) # remove kmc output suffix file

    print(f"****** Dumping & Transferring to hdfs {histFile} -> {destFile} ******")
    cmd1 = f"/usr/local/bin/kmc_dump_x {histFile} stdout | hdfs dfs -put - {destFile}"
    p = subprocess.run( cmd1, shell=True, stdout=subprocess.PIPE)

    os.remove(histFile +'.kmc_pre') # remove kmc output prefix file
    os.remove(histFile +'.kmc_suf') # remove kmc output suffix file

    return





def extractKmers( inputDataset, k, tempDir, kmcOutputPrefix):

    # run kmc on the first sequence
    # -v - verbose mode (shows all parameter settings); default: false
    # -k<len> - k-mer length (k from 1 to 256; default: 25)
    # -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12
    # -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)
    # -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9
    # -f<a/q/m/bam/kmc> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC(-fkmc); default: FASTQ
    # -ci<value> - exclude k-mers occurring less than <value> times (default: 2)
    # -cs<value> - maximal value of a counter (default: 255)
    # -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)
    # -b - turn off transformation of k-mers into canonical form
    # -r - turn on RAM-only mode
    # -n<value> - number of bins
    # -t<value> - total number of threads (default: no. of CPU cores)
    # -sf<value> - number of FASTQ reading threads
    # -sp<value> - number of splitting threads
    # -sr<value> - number of threads for 2nd stage
    # -hp - hide percentage progress (default: false)

    cmd = "/usr/local/bin/kmc -b -hp -k%d -m12 -fm -ci0 -cs1048575000 -cx2000000000 %s %s %s" % (k, inputDataset, kmcOutputPrefix, tempDir)

    print(f"****** (local) Kmer Counting {cmd} ******")

    out = subprocess.check_output(cmd.split())
    results = out.decode()
    m = re.search(r'Total no. of k-mers[ \t]*:[ \t]*(\d+)', results)
    totalKmerNumber = 0 if (m is None) else int(m.group(1))
    # print(f"****** cmd: {cmd} returned:\n{results} ******")

    m = re.search(r'No. of unique k-mers[ \t]*:[ \t]*(\d+)', results)
    totalDistinctKmerNumber = 0 if (m is None) else int(m.group(1))

    return (totalDistinctKmerNumber, totalKmerNumber)




def splitAndCount(cnt, x: str):
    cnt += 1
    return x.split('\t')[0]



def distCounter(cnt, x: str):
    cnt += 1
    return x




def countBasedMeasures(partData):
    global totKmerAAcc, totKmerBAcc, kmerStats

    D2totValue = D2ZTotValue = EuclideanTotValue = EuclideanZTotValue = 0
    Acnt = Bcnt = Ccnt = 0
    Hk1 = totalProb1 = Hk2 = totalProb2 = 0.0
    totalKmerCnt1 = float(totKmerAAcc.value)
    totalKmerCnt2 = float(totKmerBAcc.value)
    mean1 = float(kmerStats.value.mean1)
    mean2 = float(kmerStats.value.mean2)
    std1 = float(kmerStats.value.std1)
    if (std1 == 0):
        std1 = 10**-40
    std2 = float(kmerStats.value.std2)
    if (std2 == 0):
        std2 = 10**-40

    for row in partData:
        # cnt1 = 0 if row.cnt1 is None else row.cnt1
        # cnt2 = 0 if row.cnt2 is None else row.cnt2
        cnt1 = row.cnt1
        cnt2 = row.cnt2
        # zcnt1 = 0 if std1 == 0 else (cnt1 - mean1) / std1
        # zcnt2 = 0 if std2 == 0 else (cnt2 - mean2) / std2
        zcnt1 = (cnt1 - mean1) / std1
        zcnt2 = (cnt2 - mean2) / std2

        # counting based measures
        d = cnt2 - cnt1
        zd = zcnt2 -zcnt1
        EuclideanTotValue += d * d
        D2totValue +=  cnt1 * cnt2
        EuclideanZTotValue += zd * zd
        D2ZTotValue += zcnt1 * zcnt2

        # Present/Absent
        if (cnt1 > 0 and cnt2 > 0):
            Acnt += 1
        elif (cnt1 > 0):
            Bcnt += 1
        else:
            Ccnt += 1
        # calcola i valori dell'entropia in-line per evitare due pipe
        if (cnt1 > 0):
            prob1 = cnt1 / totalKmerCnt1
            totalProb1 += prob1
            Hk1 += prob1 * math.log(prob1, 2)

        if (cnt2 > 0):
            prob2 = cnt2 / totalKmerCnt2
            totalProb2 += prob2
            Hk2 += prob2 * math.log(prob2, 2)

    return iter([(EuclideanTotValue,EuclideanZTotValue, D2totValue, D2ZTotValue, Acnt, Bcnt, Ccnt, Hk1, totalProb1, Hk2, totalProb2)])




# run jaccard on sequence pair ds with kmer of length = k
def processLocalPair(seqFile1: str, seqFile2: str, k: int, theta: int, tempDir: str):
    global totDistinctKmerAAcc, totDistinctKmerBAcc, totKmerAAcc, totKmerBAcc, kmerStats

    start = time.time()

    # first locally extract kmer statistics for both sequences

    baseSeq1 = Path(seqFile1).stem
    kmcOutputPrefixA = f"{tempDir}/{baseSeq1}-k={k}"
    destFilenameA = f"{hdfsDataDir}/{baseSeq1}-k={k}.txt"
    # calcola comunque kmc per avere i valori di totDistinctKmerA, totKmerA
    (totDistinctKmerA, totKmerA) = extractKmers(seqFile1, k, tempDir, kmcOutputPrefixA)
    if (not checkPathExists(destFilenameA)):
        # load kmers statistics from histogram files (dumping kmc output to hdfs)
        loadHistogramOnHDFS(kmcOutputPrefixA, destFilenameA)
    else:
        # altrimenti rimuove solo i file temporanei di KMC e usera' destFilenameA come input
        os.remove(kmcOutputPrefixA+'.kmc_pre')
        os.remove(kmcOutputPrefixA+'.kmc_suf')

    baseSeq2 = Path(seqFile2).stem
    kmcOutputPrefixB = f"{tempDir}/{baseSeq2}-k={k}"
    destFilenameB = f"{hdfsDataDir}/{baseSeq2}-k={k}.txt"
    # calcola comunque kmc per avere i valori di totDistinctKmerB, totKmerB
    (totDistinctKmerB, totKmerB) = extractKmers(seqFile2, k, tempDir, kmcOutputPrefixB)
    if (not checkPathExists(destFilenameB)):
        # load kmers statistics from histogram files (dumping kmc output to hdfs)
        loadHistogramOnHDFS(kmcOutputPrefixB, destFilenameB)
    else:
        # altrimenti rimuove solo i file temporanei di KMC e usera' destFilenameB come input
        os.remove(kmcOutputPrefixB+'.kmc_pre')
        os.remove(kmcOutputPrefixB+'.kmc_suf')

    #
    # inizio procedura Dataframe oriented (out of memory)
    #
    # Broadcast variables ... for all map tasks. These are read-only
    totDistinctKmerAAcc = sc.broadcast(totDistinctKmerA)
    totDistinctKmerBAcc = sc.broadcast(totDistinctKmerB)
    totKmerAAcc = sc.broadcast(totKmerA)
    totKmerBAcc = sc.broadcast(totKmerB)

    schema1 = StructType([ StructField('kmer', StringType(), True), StructField('cnt1', IntegerType(), True)])
    schema2 = StructType([ StructField('kmer', StringType(), True), StructField('cnt2', IntegerType(), True)])

    df1 = spark.read.format("csv").schema(schema1).options(delimiter='\t').load(destFilenameA)
    df2 = spark.read.format("csv").schema(schema2).options(delimiter='\t').load(destFilenameB)

    # inner solo l'intersezione
    # inner = df1.join(df2, df1.kmer == df2.kmer, "inner")
    # full outer join di tutte le righe (poiche' le colonne chiave hanno lo stesso nome possiamo usare una lista
    # al termine tutti i valori null (i valori di cnt per k-mer non definiti) vengono posti a 0
    outer = df1.join(df2, ['kmer'], how='full').na.fill(value=0, subset=['cnt1','cnt2'])
    # questa produce 2 colonne kmer
    # outer = df1.join(df2, df1.kmer == df2.kmer, "outer")

    # df4 = outer.select(EuclidUDF(col('cnt1'), col('cnt2'))).show()

    # calcola media e deviazione standard delle frequenze.
    stats = outer.select(sf.mean(sf.col('cnt1')).alias('mean1'),
                         sf.stddev(sf.col('cnt1')).alias('std1'),
                         sf.mean(sf.col('cnt2')).alias('mean2'),
                         sf.stddev(sf.col('cnt2')).alias('std2')
                         ).collect()

    # share the results to all workers
    kmerStats = sc.broadcast( stats[0])
    allDist = outer.rdd.mapPartitions(countBasedMeasures).collect()

    totEuclidZ = totD2Z = 0.0
    totEuclid = totD2 = Acnt = Bcnt = Ccnt = 0
    HkA = totalProbA = HkB = totalProbB = 0.0
    for t in allDist:
        totEuclid += t[0]; totEuclidZ += t[1]; totD2 += t[2]; totD2Z += t[3]
        Acnt += t[4]; Bcnt += t[5]; Ccnt += t[6]
        HkA  += t[7]; totalProbA += t[8]; HkB  += t[9]; totalProbB += t[10]

    vec = [0,0,0,0,0,0,0,0,0,0,0]
    for t in allDist:
        vec = [ x[0] + x[1] for x in zip(vec, t) ]

    print(f"****** {vec[0]}, {vec[1]}, {vec[2]}, {vec[3]}, {vec[4]}, {vec[5]}, {vec[6]}, {vec[7]}, {vec[8]} ******")

    if round(totalProbA,0) != 1.0:
        # raise ValueError("Somma(Pa = {round(totalProbA, 0):f} must be 1.0. Aborting")
        print(f"****** Somma(Pa) = {round(totalProbA, 0):.2f} must be 1.0!!! ******")

    entropySeqA = EntropyData( totDistinctKmerA, totKmerA, HkA)

    if round(totalProbB,0) != 1.0:
        # raise ValueError("Somma(Pb) = {round(totalProbB, 0):f} must be 1.0. Aborting")
        print(f"****** Somma(Pb) = {round(totalProbB, 0):.2f} must be 1.0!!! ******")

    entropySeqB = EntropyData( totDistinctKmerB, totKmerB, HkB)

    euclideanDistance = math.sqrt(totEuclid)
    euclideanDistanceZ = math.sqrt(totEuclidZ)
    print(f"****** Euclidean = {euclideanDistance:.4f}, EuclideanZ = {euclideanDistanceZ:.4f} ******")
    print(f"****** D2 = {totD2:,} D2Z = {totD2Z:.4f} ******")
    print(f"****** Present/Absent = {Acnt:,}, {Bcnt:,}, {Ccnt:,} ******")
    print(f"****** HkA: {entropySeqA.Hk:.5f} HkB: {entropySeqB.Hk:.5f}, totDistinctKmerA: {entropySeqA.totalKmerCnt:,}, totDistinctKmerB: {entropySeqB.totalKmerCnt:,} ******")

# dati3 = runCountBasedMeasures(cnts, k)
    dati3 =  [totD2, totD2Z, euclideanDistance, euclideanDistanceZ]

    # implementazione precedente per avere A, B, e C. Effettuato test i risultati coincidono
    # tot1Acc = sc.accumulator(0)
    # seq1 = sc.textFile(destFilenameA).map( lambda x: splitAndCount( tot1Acc, x))
    #
    # tot2Acc = sc.accumulator(0)
    # seq2 = sc.textFile(destFilenameB).map( lambda x: splitAndCount( tot2Acc, x))
    #
    # intersection = seq1.intersection(seq2)
    #
    # outputFile = '%s/k=%d-%s-%s.txt' % (hdfsDataDir, k, baseSeq1, baseSeq2)
    #
    # # tot3Acc = sc.accumulator(0)
    # # intersection.map(lambda x: distCounter(tot3Acc, x)).saveAsTextFile( outputFile)
    # # bothCnt = tot3Acc.value
    #
    # ACnt = intersection.count()
    # Bcnt = tot1Acc.value - bothCnt
    # Ccnt = tot2Acc.value - bothCnt

    ###################################################
    dati1 = runPresentAbsent(Acnt, Bcnt, Ccnt, k)

    dati2 = runMash(seqFile1, seqFile2, k)

    # dati4 = dati entropia
    dati4 = entropySeqA.toString() + entropySeqB.toString()

    delay = time.time()-start
    dati0 = [baseSeq1, baseSeq2, start, delay, theta, k]

    # do not remove base textual histogram file (for all k values) because can be reused in next iterations (for other theta values)
    # cmd = f"hdfs dfs -rm -skipTrash {destFilenameA}"
    # p = subprocess.Popen(cmd.split())
    # p.wait()

    # remove textual histogram files from hdfs
    # cmd = f"hdfs dfs -rm -skipTrash {destFilenameB}"
    # p = subprocess.Popen(cmd.split())
    # p.wait()
    
    return dati0 + dati1 + dati2 + dati3 + dati4    # nuovo record output




def writeHeader( writer):#
    columns0 = ['sequenceA', 'sequenceB', 'start time', 'real time', 'Theta', 'k'] # dati 0
    columns1 = [ 'A', 'B', 'C', 'D', 'N', 'A/N',
               'Anderberg', 'Antidice', 'Dice', 'Gower', 'Hamman', 'Hamming',
                'Jaccard', 'Kulczynski', 'Matching', 'Ochiai',
                'Phi', 'Russel', 'Sneath', 'Tanimoto', 'Yule']

    columns2 = []
    for ss in sketchSizes:
        columns2.append( 'Mash Pv (%d)' % ss)
        columns2.append( 'Mash Distance(%d)' % ss)
        columns2.append( 'A (%d)' % ss)
        columns2.append( 'N (%d)' % ss)

    columns3 = [ 'D2', 'D2Z', 'Euclidean', 'EuclideanZ']

    columns4 = ['NKeysA', '2*totalCntA', 'deltaA', 'HkA', 'errorA',
                'NKeysB', '2*totalCntB', 'deltaB', 'HkB', 'errorB']

    writer.writerow(columns0 + columns1 + columns2 + columns3 + columns4)
    
                 




# processa localmente una coppia di sequenze seqFile1 e seqFile2
def processPairs(seqFile1: str, seqFile2: str, theta: int):
    # process local sequence files in the same local directory (temporary named ttt)
    tempDir = os.path.dirname( seqFile1)+'/ttt'
    if (not os.path.isdir(tempDir)):
        os.mkdir(tempDir)

    # local file system result file
    outFile = "%s/%s-%s-%02d-%d.csv" % (os.path.dirname( seqFile1), Path(seqFile1).stem,Path(seqFile2).stem, theta, int(time.time()))
    with open(outFile, 'w') as file:
        csvWriter = csv.writer(file)        
        writeHeader(csvWriter)
        file.flush()

        if (seqFile2 == "synthetic"):
            # produce il file allontanato da seqFile1 di un fattore theta
            (f, ext) = os.path.splitext(seqFile1)
            seqFile2 = f"{f}-{theta}{ext}"
            mkd.MoveAwaySequence(seqFile1, seqFile2, theta)

        for k in range( minK, maxK+1, stepK):
            # run kmc on both the sequences and eval A, B, C, D + Mash + Entropy
            print(f"****** Starting {Path(seqFile1).stem} vs {Path(seqFile2).stem} k = {k} T = {theta} ******")
            res = processLocalPair(seqFile1, seqFile2, k, theta, tempDir)
            csvWriter.writerow( res)
            file.flush()

            
    # clean up
    # do not remove dataset on hdfs
    # remove histogram files (A & B) + mash sketch file and kmc temporary files
    try:
        print(f"****** Cleaning temporary directory {tempDir} ******")
        shutil.rmtree(tempDir)
    except OSError as e:
        print("Error removing: %s: %s" % (tempDir, e.strerror))






def main():
    global hdfsDataDir, hdfsPrefixPath, spark, sc, thetaValue, minK, maxK


    hdfsDataDir = hdfsPrefixPath

    argNum = len(sys.argv)
    if (argNum < 5 or argNum > 6):
        """
            Usage: PySparkPASingleSequenceOutMemory Sequence1 Sequence2 theta dataDir [kValue]
        """
    else:
        # theta viene utilizzato SOLO se Sequence2 == "synthetic" altrimenti viene ignorato
        thetaValue = int(sys.argv[3])
        hdfsDataDir = '%s/%s' % (hdfsPrefixPath, sys.argv[4])

    if (argNum == 6):
        # use just one k value instead of all values from minK to maxK (step)
        minK = int(sys.argv[5])
        maxK = int(sys.argv[5])

    seqFile1 = sys.argv[1] # le sequenze sono sul file system locale
    seqFile2 = sys.argv[2] # per eseguire localmente l'estrazione dei k-mers
    # outFile = '%s/%s-%s.csv' % (hdfsDataDir, Path( seqFile1).stem, Path(seqFile2).stem )

    if (seqFile2 == "synthetic"):
        print(f"****** Comparing: {Path(seqFile1).stem} vs {Path(seqFile2).stem} with {minK} <= k <= {maxK} and Theta = {thetaValue} in hdfsDataDir = {hdfsDataDir} ******")
    else:
        print(f"****** Comparing: {Path(seqFile1).stem} vs {Path(seqFile2).stem} with {minK} <= k <= {maxK} in hdfsDataDir = {hdfsDataDir} ******")

    spark = SparkSession \
        .builder \
        .appName( f"{Path( sys.argv[0]).stem} {Path(seqFile1).stem} {Path(seqFile2).stem} {minK} <= k <= {maxK} theta = {thetaValue}") \
        .getOrCreate()

    sc = spark.sparkContext

    sc2 = spark._jsc.sc()
    nWorkers =  len([executor.host() for executor in sc2.statusTracker().getExecutorInfos()]) - 1

    if (not checkPathExists( hdfsDataDir)):
        print(f"****** Data dir: {hdfsDataDir} does not exist. Program terminated. ******")
        exit( -1)

    print(f"****** {nWorkers} workers, hdfsDataDir: {hdfsDataDir} ******")

    processPairs(seqFile1, seqFile2, thetaValue)

    spark.stop()


    



if __name__ == "__main__":
    main()
