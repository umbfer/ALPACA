#! /usr/bin/python3
import numpy
import re
import os
import sys
import glob
import tempfile
import shutil
import subprocess
import csv
import math
from filelock import Timeout, FileLock
import numpy as np

from operator import add
import pyspark
from pyspark.sql import SparkSession
from pyspark import SparkFiles


hdfsPrefixPath = 'hdfs://master2:9000/user/cattaneo/data/'
spark = None

models = ['Uniform', 'MotifRepl-U', 'PatTransf-U', 'Uniform-T1']
hdfsDataDir = hdfsPrefixPath + 'dataset5-1000'
inputRE = '*1000.[2-8]00000*'
#lengths = range(1000, 50001, 1000) # small dataset
lengths = [ 10000, 100000, 1000000, 10000000]
gVals = [10, 50, 100]
nTests = 500
minK = 4
maxK = 33
stepK = 4
sketchSize = 1000
outFile = 'PresentAbsentData.csv'






# calcola anche i valori dell'entropia per non caricare due volte l'istogramma
def loadKmerList( file):

    print("Loading from file %s" % file)
    totalCnt = 0
    # ogni file contiene l'istogramma di una sola sequenza prodotto con kmc 3
    with open(file) as inFile:
        seqDict = dict()
        for line in inFile:
            s = line.split()   # molto piu' veloce della re
            if (len(s) != 2):
                print( "%sMalformed histogram file (%d token)" % (line, len(s)))
                exit()
            else:
                kmer = s[0]
                count = int(s[1])

            if kmer in seqDict:
                seqDict[kmer] = seqDict[kmer] + count
            else:
                seqDict[kmer] = count

            totalCnt = totalCnt + count

    totalKmers = 0 # this should be seqLen - k + 1
    totalKeys = 0  # this value should be len(seqDict.keys()) 
    totalProb = 0.0
    Hk = 0.0
    kList = seqDict.keys()
    aSize = len(kList)
    npArray = np.empty( aSize, numpy.str)
    for key in kList:
        cnt = seqDict[key]
        prob = cnt / float(totalCnt)
        totalProb = totalProb + prob
        npArray[totalKeys] = key
        totalKeys = totalKeys + 1
        totalKmers = totalKmers + cnt
        Hk = Hk + prob * math.log(prob, 2)
        # print( "prob(%s) = %f log(prob) = %f" % (key, prob, math.log(prob, 2)))

    Hk = Hk * -1
    # arr = np.array(list(seqDict.keys()))
    nKeys = npArray.size
    if (totalKeys != nKeys):
        raise ValueError( "TotalKeys = %d vs len(seqDict.keys()) = %d" % (totalKmers, nKeys))

    if (totalKmers != totalCnt):
        raise ValueError ( "k-mers count %d vs %d (should be = seqLen - k + 1)" % (totalKmers, totalCnt))

    if (round(totalProb,0) != 1.0):
        raise ValueError("Somma(p) = %f must be 1.0. Aborting" % round(totalProb, 0))

    return (nKeys, totalCnt, Hk, npArray)



def extractKmers( dataset, k, seq, seqDir, tempDir):

    inputDataset = '%s/%s-%s.fasta' % (seqDir, dataset, seq)
    kmcOutputPrefix = "%s/k=%d%s-%s" % (tempDir, k, dataset, seq)
    # run kmc on the first sequence
    cmd = "/usr/local/bin/kmc -b -k%d -m2 -fm -ci0 -cs1000000 %s %s %s" % (k, inputDataset, kmcOutputPrefix, tempDir)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))

    # dump the result -> kmer histogram
    histFile = "%s/distk=%d_%s-%s.hist" % (tempDir, k, dataset, seq)
    cmd = "/usr/local/bin/kmc_dump %s %s" % ( kmcOutputPrefix, histFile)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))
    # load kmers from histogram file
    results = loadKmerList(histFile)

    return results




# run jaccard on sequence pair ds with kmer of length = k
def runPresentAbsent( ds, model, seqId, seqLen, gamma, k, seqDir, tempDir):

    m = re.search(r'^(.*)-(\d+)\.(\d+)(.*)', ds)
    if (m is not None):
        if ((model != m.group(1)) or (seqId != int(m.group(2))) or
            (seqLen != int(m.group(3)))):
            # or (gamma != float('0.' + m.group(4)[5:]))):
            err = 'Wrong function paramenters: %s %s, %d %d, %d %d, %s %s' % (
                model, m.group(1), seqId, int(m.group(2)),
                seqLen, int(m.group(3)), gamma, float('0.' + m.group(4)[5:]))
            return [ err]
    else:
        err = 'Malformed dataset name'
        return [ err]

    (nKeys, totalCnt, Hk, leftKmers) = extractKmers(ds, k, 'A', seqDir, tempDir)
    (nKeys, totalCnt, Hk, rightKmers) = extractKmers(ds, k, 'B', seqDir, tempDir)
    print("left: %d, right: %d" % (leftKmers.size, rightKmers.size))

    intersection = np.intersect1d( leftKmers, rightKmers)
    bothCnt = intersection.size
    A = bothCnt
    leftCnt = leftKmers.size - bothCnt
    B = leftCnt
    rightCnt = rightKmers.size - bothCnt
    C = rightCnt
    
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
        anderberg = 'UnDefined'

    # Antidice dissimilarity => Antidice = 1 - A/(A + 2(B + C))
    try:
        antidice = 1 - A / float(A + 2.0 * (B + C))
    except (ZeroDivisionError, ValueError):
        antidice = 'UnDefined'

    # Dice dissimilarity => Dice = 1 - 2A/(2A + B + C)
    try:
        dice = 1 - 2*A / float(2.0*A + B + C)
    except (ZeroDivisionError, ValueError):
        dice  = 'UnDefined'
    # Gower dissimilarity => Gower = 1 - A x D/sqrt(A + B) x(A + C) x (D + B x (D + C)
    try:
        gower = 1 - A * D / math.sqrt((A + B) * (A + C) * (D + B * (D + C)))
    except (ZeroDivisionError, ValueError):
        gower = 'UnDefined'

    # Hamman dissimilarity => Hamman = 1 - [((A + D) - (B + C))/N]2
    try:
        hamman = 1 - math.pow((((A + D) - (B + C)) / float(NMax)), 2.0)
    except (ZeroDivisionError, ValueError):
        hamman = 'UnDefined'

    # Hamming dissimilarity => Hamming = (B + C)/N
    try:
        hamming = (B + C)/ NMax
    except (ZeroDivisionError, ValueError):
        hamming = 'UnDefined'

    # Jaccard dissimilarity => Jaccard = 1 - A/(N - D)
    try:
        jaccard = 1 - A / (NMax - D)
    except (ZeroDivisionError, ValueError):
        jaccard = 'UnDefined'

    jaccardDistance = 1 - min( 1.0, A / float(NMax - D))

    # Kulczynski dissimilarity => Kulczynski = 1 - (A/(A + B) + A/(A + C)) / 2
    try:
        kulczynski = 1 - (A / float(A + B) + A / float(A + C)) / 2.0
    except (ZeroDivisionError, ValueError):
        kulczynski = 'UnDefined'

    # Matching dissimilarity => Matching = 1 - (A + D)/N
    try:
        matching = 1 - (A + D) / NMax
    except (ZeroDivisionError,ValueError):
        matching = 'UnDefined'

    # Ochiai dissimilarity => Ochiai = 1 - A/sqrt(A + B) x (A + C)
    try:
        ochiai = 1 - A / math.sqrt((A + B) * (A + C))
    except (ZeroDivisionError, ValueError):
        ochiai = 'UnDefined'

    # Phi dissimilarity => Phi = 1 - [(A x  B x  C x D)/sqrt(A + B) x (A + C) x (D + B) x (D + C)]2
    try:
        phi = 1 - math.pow((A * B * C * D)/ math.sqrt((A + B) * (A + C) * (D + B) * (D + C)), 2.0)
    except (ZeroDivisionError, ValueError):
        phi = 'UnDefined'

    # Russel dissimilarity => Russel = 1 - A/N
    try:
        russel = 1 - A / NMax
    except (ZeroDivisionError, ValueError):
        russel = 'UnDefined'

    # Sneath dissimilarity => Sneath = 1 - 2(A + D)/(2 x (A + D) + (B + C))
    try:
        sneath = 1 - 2.0 * (A + D) / (2.0 * (A + D) + (B + C))
    except (ZeroDivisionError, ValueError):
        sneath = 'UnDefined'

    # Tanimoto dissimilarity => Tanimoto = 1 - (A + D)/((A + D) + 2(B + C))
    try:
        tanimoto = 1 - (A + D) / float((A + D) + 2.0 * (B + C))
    except (ZeroDivisionError, ValueError):
        tanimoto = 'UnDefined'

    # Yule dissimilarity => Yule = 1 - [(A x D - B x C)/(A x D + B x C)]2
    try:
        yule = 1 - math.pow(((A * D - B * C) / float(A * D + B * C)), 2.0)
    except (ZeroDivisionError, ValueError):
        yule = 'UnDefined'

    # run mash on the same sequence pair
    inputDS1 = '%s/%s-A.fasta' % (seqDir, ds)
    inputDS2 = '%s/%s-B.fasta' % (seqDir, ds)

    # extract mash sketch from the first sequence
    cmd = "/usr/local/bin/mash sketch -s %d -k %d %s" % (sketchSize, k, inputDS1)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))

    cmd = "/usr/local/bin/mash sketch -s %d -k %d %s" % (sketchSize, k, inputDS2)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))

    cmd = "/usr/local/bin/mash dist %s.msh %s.msh" % (inputDS1, inputDS2)
    out = subprocess.check_output(cmd.split())
    print("cmd: %s returned: %s" % (cmd, out))
    mashResults = out.split()
    mashPv = float(mashResults[2])
    mashDist = float(mashResults[3])
    mashAN = mashResults[4].decode('UTF-8')
    # dati sull'entropia della sequenza B (non A)
    delta = float(nKeys) / (2 * totalCnt)

    # salva il risultato nel file CSV
    data = [model, gamma, seqLen, seqId, k, A, B, C, D, NMax,
            anderberg, antidice, dice, gower, hamman, hamming, jaccard, jaccardDistance,
            kulczynski, matching, ochiai, phi, russel, sneath, tanimoto, yule,
            mashPv, mashDist, mashAN,
            nKeys, 2 * totalCnt, delta, Hk, delta / Hk]

    # clean up remove kmc temporary files
    for f in glob.glob('%s/*%s-*' % (tempDir, ds)):
           os.remove(f)
    os.remove(inputDS1 + '.msh')
    os.remove(inputDS2 + '.msh')

    return data





def splitFastaSequences(model, seqLen, gamma, dsContent, nSeq):

    seqDir = tempfile.mkdtemp()

    lookForHeader = True
    start = 0
    cnt = 0
    for ndx in range(len( dsContent)):
        if (dsContent[ndx] == 10):  # end of line
            if lookForHeader:
                header = dsContent[start:ndx+1].decode('UTF-8') #incluso il '\n'
                start = ndx + 1
                lookForHeader = False
                m = re.search(r'^>(.+)\.(\d+)(.*)-([AB]$)', header)
                if (m is not None):
                    # seqName = m.group(1) # unico per l'intero file
                    seqId = int(m.group(2))
                    # gValue = m.group(3) # unico per l'intero file
                    pairId =  m.group(4)
                    # print( "Name: %s, id:%d, GValue: %s, pair:%s" % (model, seqId, gamma, pairId))
                else:
                    print("Malformed header: %s" % header)
                    return ''
            else:
                sequence = dsContent[start:ndx+1].decode('UTF-8') # incluso il '\n'
                start = ndx + 1
                lookForHeader = True
                # salva la sequenza localmente
                fileName = "%s/%s-%04d.%d%s-%s.fasta" % (seqDir, model, seqId, seqLen, gamma, pairId)
                with open(fileName, "w") as outText:
                    outText.write("%s%s" % (header, sequence)) # \n are in the original strings

                cnt = cnt + 1
                if (cnt > nSeq*2):    # salva solo le prime nSeq coppie che servono
                    break


    return seqDir




def processDataset( ds):

    m = re.search(r'^(.*)-(\d+)\.(\d+)(.*).fasta', os.path.basename(ds[0]))
    if (m is None):
        print("Malformed file name %s" % ds[0])
        return -2
    else:
        model = m.group(1)
        nPairs = int(m.group(2))
        seqLen = int(m.group(3))
        gamma = m.group(4)

    if (nTests > nPairs):
        raise ValueError( 'For this dataset the max number of runs is %d.' % nPairs)

    print("Splitting dataset: %s@%s" % (ds[0], os.uname()[1]))
    seqDir = splitFastaSequences( model, seqLen, gamma, ds[1], nTests)

    # process private temporary directory
    tempDir = tempfile.mkdtemp()

    results = []
    for seqId in range(1, nTests+1):
        ds = '%s-%04d.%d%s' % (model, seqId, seqLen, gamma)
        for k in range( minK, maxK, stepK):
            # run kmc on both the sequences and eval A, B, C, D + Mash + Entropy
            results.append(runPresentAbsent(ds, model, seqId, seqLen, gamma, k, seqDir, tempDir))

    # clean up
    # do not remove daset on hdfs
    # remove histogram files (A & B) + mash sketch file and kmc temporary files
    try:
        print("Cleaning directory %s and %s" % (tempDir, seqDir))
        shutil.rmtree(tempDir)
        shutil.rmtree(seqDir)
    except OSError as e:
        print("Error removing: %s and %s: %s" % (tempDir, seqDir, e.strerror))

    return results



def main():
    global spark
    """
       Usage: PySparkPresentAbsent [partitions]
    """

    spark = SparkSession \
        .builder \
        .appName("PresentAbsent") \
        .getOrCreate()

    sc = spark.sparkContext

    # partitions = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    #
    rdd = sc.binaryFiles( '%s/%s.fasta' % (hdfsDataDir, inputRE))
    # print("Number of Partitions: " + str(rdd.getNumPartitions()))
    #
    counts = rdd.flatMap(lambda x: processDataset(x))

    columns = ['model', 'gamma', 'seqLen', 'pairId', 'k', 'A', 'B', 'C', 'D', 'N',
                'Anderberg', 'Antidice', 'Dice', 'Gower', 'Hamman', 'Hamming',
                'Jaccard', 'jaccardDistance', 'Kulczynski', 'Matching', 'Ochiai',
                'Phi', 'Russel', 'Sneath', 'Tanimoto', 'Yule',
                'Mash Pv', 'Mash Distance', 'A/N',
                'NKeys', '2*totalCnt', 'delta', 'Hk', 'error']
    df = counts.toDF(columns)
    # df.show()
    # print( counts.collect())

    # df.write.format("csv").save(outFile)
    df.write.option("header",True).csv(outFile)
    spark.stop()






if __name__ == "__main__":
    main()
