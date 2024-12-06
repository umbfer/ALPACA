#! /usr/bin/python

import re
import os
import sys
import glob
import subprocess
import csv
import math
from filelock import Timeout, FileLock
import numpy as np

libPath = '/home/cattaneo/spark/power_statistics/scripts'
if (os.path.exists(libPath)):
    sys.path.append(libPath)
    import splitFasta

    
# process private temporary directory
tempDir = "tmp.%d" % os.getpid()

models = ['Uniform', 'MotifRepl-U', 'PatTransf-U', 'Uniform-T1']
hdfsDataDir = 'data/dataset5-1000'
#lengths = range(1000, 50001, 1000) # small dataset
lengths = [ 1000000, 10000000]
gVals = [10, 50, 100]
nPairs = 1000
nTests = 500
minK = 4
maxK = 32
sketchSize = 1000
outFile = 'PresentAbsentData.csv'

# variabili globali per il calcolo dell'entropia
Hk = 0.0
nKeys = 0
totalCnt = 0


def extractKmers( dataset, k, seq):
    inputDataset = '%s/%s-%s.fasta' % (splitFasta.seqDistDir, dataset, seq)
    kmcOutputPrefix = "%s/k=%d%s-%s" % (tempDir, k, dataset, seq)
    # run kmc on the first sequence
    cmd = "kmc -b -k%d -m2 -fm -ci0 -cs1000000 %s %s %s" % (k, inputDataset, kmcOutputPrefix, tempDir)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))

    # dump the result -> kmer histogram
    histFile = "%s/distk=%d_%s-%s.hist" % (tempDir, k, dataset,seq)
    cmd = "kmc_dump %s %s" % ( kmcOutputPrefix, histFile)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))
    # load kmers from histogram file
    vect = loadKmerList(histFile)
    # remove temporary files
    # os.remove(histFile)
    # for f in glob.glob(kmcOutputPrefix+'*'):
    #    os.remove(f)
    return vect




# calcola anche i valori dell'entropia per non caricare due volte l'istogramma
def loadKmerList( file):
    global Hk, nKeys, totalCnt
    
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
    for key in seqDict.keys():
        cnt = seqDict[key]
        prob = cnt / float(totalCnt)
        totalProb = totalProb + prob
        totalKeys = totalKeys + 1
        totalKmers = totalKmers + cnt
        Hk = Hk + prob * math.log(prob, 2)
        # print( "prob(%s) = %f log(prob) = %f" % (key, prob, math.log(prob, 2)))

    Hk = Hk * -1
    nKeys = len(seqDict.keys())
    if (totalKeys != nKeys):
        print( "errore totalKeys = %d vs len(seqDict.keys()) = %d" % (totalKmers, nKeys))
        exit(-1)


    if (totalKmers != totalCnt):
        print( "errore k-mers count %d vs %d (this should be = seqLen - k + 1)" % (totalKmers, totalCnt))
        exit(-1)

    if (round(totalProb,0) != 1.0):
        print( "errore Somma(p) = %f" % (totalProb))
        exit(-1)

    return np.array(seqDict.keys())



# run jaccard on sequence pair ds with kmer of length = k
def runPresentAbsent( ds, k):
    global Hk, nKeys, totalCnt, tempDir

    m = re.search(r'^(.*)-(\d+)\.(\d+)(.*)', ds)
    if (m is not None):
        model = m.group(1)
        pairId = int(m.group(2))
        seqLen = int(m.group(3))
        gamma = float('0.' + m.group(4)[5:])
    else:
        print('Malformed dataset name')


    leftKmers = extractKmers(ds, k, 'A')
    rightKmers = extractKmers(ds, k, 'B')

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
    except ZeroDivisionError:
        anderberg = 'UnDefined'

    # Antidice dissimilarity => Antidice = 1 - A/(A + 2(B + C))
    try:
        antidice = 1 - A / float(A + 2.0 * (B + C))
    except ZeroDivisionError:
        antidice = 'UnDefined'

    # Dice dissimilarity => Dice = 1 - 2A/(2A + B + C)
    try:
        dice = 1 - 2*A / float(2.0*A + B + C)
    except ZeroDivisionError:
        dice  = 'UnDefined'
    # Gower dissimilarity => Gower = 1 - A x D/sqrt(A + B) x(A + C) x (D + B x (D + C)
    try:
        gower = 1 - A * D / math.sqrt((A + B) * (A + C) * (D + B * (D + C)))
    except ZeroDivisionError:
        gower = 'UnDefined'

    # Hamman dissimilarity => Hamman = 1 - [((A + D) - (B + C))/N]2
    try:
        hamman = 1 - math.pow((((A + D) - (B + C)) / float(NMax)), 2.0)
    except ZeroDivisionError:
        hamman = 'UnDefined'

    # Hamming dissimilarity => Hamming = (B + C)/N
    try:
        hamming = (B + C)/ NMax
    except ZeroDivisionError:
        hamming = 'UnDefined'

    # Jaccard dissimilarity => Jaccard = 1 - A/(N - D)
    try:
        jaccard = 1 - A / (NMax - D)
    except ZeroDivisionError:
        jaccard = 'UnDefined'

    jaccardDistance = 1 - min( 1.0, A / float(NMax - D))

    # Kulczynski dissimilarity => Kulczynski = 1 - (A/(A + B) + A/(A + C)) / 2
    try:
        kulczynski = 1 - (A / float(A + B) + A / float(A + C)) / 2.0
    except ZeroDivisionError:
        kulczynski = 'UnDefined'

    # Matching dissimilarity => Matching = 1 - (A + D)/N
    try:
        matching = 1 - (A + D) / NMax
    except ZeroDivisionError:
        matching = 'UnDefined'

    # Ochiai dissimilarity => Ochiai = 1 - A/sqrt(A + B) x (A + C)
    try:
        ochiai = 1 - A / math.sqrt((A + B) * (A + C))
    except ZeroDivisionError:
        ochiai = 'UnDefined'

    # Phi dissimilarity => Phi = 1 - [(A x  B x  C x D)/sqrt(A + B) x (A + C) x (D + B) x (D + C)]2
    try:
        phi = 1 - math.pow((A * B * C * D)/ math.sqrt((A + B) * (A + C) * (D + B) * (D + C)), 2.0)
    except ZeroDivisionError:
        phi = 'UnDefined'

    # Russel dissimilarity => Russel = 1 - A/N
    try:
        russel = 1 - A / NMax
    except ZeroDivisionError:
        russel = 'UnDefined'

    # Sneath dissimilarity => Sneath = 1 - 2(A + D)/(2 x (A + D) + (B + C))
    try:
        sneath = 1 - 2.0 * (A + D) / (2.0 * (A + D) + (B + C))
    except ZeroDivisionError:
        sneath = 'UnDefined'

    # Tanimoto dissimilarity => Tanimoto = 1 - (A + D)/((A + D) + 2(B + C))
    try:
        tanimoto = 1 - (A + D) / float((A + D) + 2.0 * (B + C))
    except ZeroDivisionError:
        tanimoto = 'UnDefined'

    # Yule dissimilarity => Yule = 1 - [(A x D - B x C)/(A x D + B x C)]2
    try:
        yule = 1 - math.pow(((A * D - B * C) / float(A * D + B * C)), 2.0)
    except ZeroDivisionError:
        yule = 'UnDefined'

    # run mash on the same sequence pair
    inputDS1 = '%s/%s-A.fasta' % (splitFasta.seqDistDir, ds)
    inputDS2 = '%s/%s-B.fasta' % (splitFasta.seqDistDir, ds)

    # run kmc on the first sequence
    cmd = "mash sketch -s %d -k %d %s" % (sketchSize, k, inputDS1)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))
    cmd = "mash sketch -s %d -k %d %s" % (sketchSize, k, inputDS2)
    p = subprocess.Popen(cmd.split())
    p.wait()
    print("cmd: %s returned: %s" % (cmd, p.returncode))
    cmd = "mash dist %s.msh %s.msh" % (inputDS1, inputDS2)
    out = subprocess.check_output(cmd.split())
    print("cmd: %s returned: %s" % (cmd, out))

    mashResults = out.split()

    # dati sull'entropia della sequenza B (non A)
    delta = float(nKeys) / (2 * totalCnt)

    # salva il risultato nel file CSV
    data = [model, gamma, seqLen, pairId, k, A, B, C, D, NMax,
            anderberg, antidice, dice, gower, hamman, hamming, jaccard, jaccardDistance,
            kulczynski, matching, ochiai, phi, russel, sneath, tanimoto, yule,
            mashResults[2], mashResults[3], mashResults[4],
            nKeys, 2*totalCnt, delta, Hk, delta/Hk ]

    lock = FileLock(outFile + '.lck')
    try:
        lock.acquire(5)
        print("Lock acquired.")
        print( data)
        if (not os.path.exists( outFile)):
            f = open(outFile, 'w')
            writer = csv.writer(f)
            header = ['model', 'gamma', 'seqLen', 'pairId', 'k', 'A', 'B', 'C', 'D', 'N',
                      'Anderberg', 'Antidice', 'Dice', 'Gower', 'Hamman', 'Hamming',
                      'Jaccard', 'jaccardDistance', 'Kulczynski', 'Matching', 'Ochiai',
                      'Phi', 'Russel', 'Sneath', 'Tanimoto', 'Yule',
                      'Mash Pv', 'Mash Distance', 'A/N',
                      'NKeys', '2*totalCnt', 'delta', 'Hk', 'error']
            writer.writerow(header)
        else:
            f = open(outFile, 'a')
            writer = csv.writer(f)

        writer.writerow(data)
        f.close()
        lock.release()

    except Timeout:
        print("Another instance of this application currently holds the lock.")

    # clean up remove kmc temporary files
    for f in glob.glob('%s/*%s-*' % (tempDir, ds)):
            os.remove(f)



def main():

    # process private temporary directory
    os.mkdir(tempDir)

    for seqLen in lengths:
        for model in models:
            gammas = [ ' ' ]  if model.startswith('Uniform') else gVals        
            for g in gammas:
                # PatTransf-U-1000.9600000.G=0.100.fasta oppure Uniform-1000.1000000.fasta
                gamma = '' if g == ' ' else '.G=0.%03d' % g
                dataset = '%s-%d.%d%s.fasta' % (model, nPairs, seqLen, gamma)
                
                # download from hdfs the dataset
                cmd = "hdfs dfs -get %s/%s ." % (hdfsDataDir, dataset) 
                p = subprocess.Popen(cmd.split())
                p.wait()
                print("cmd: %s returned: %s" % (cmd, p.returncode))
                
                splitFasta.splitFastaSequences( dataset)
            
                for seqId in range(1, nTests+1):
                    
                    ds = '%s-%04d.%d%s' % (model, seqId, seqLen, gamma)
                    for k in range(4, 33, 4):
                        
                        # run kmc on both the sequences and eval A, B, C, D + Mash + Entropy
                        runPresentAbsent(ds, k)

                # clean up
                # posso rimnuovere il dataset se importato altrimenti non cancellare !!!!
                os.remove( dataset)
                # remove histogram files (A & B) + mash sketch files
                for f in glob.glob('%s/%s-*' % (splitFasta.seqDistDir, model)):
                    os.remove(f)




                        

if __name__ == "__main__":
    main()
