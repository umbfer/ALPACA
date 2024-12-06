#! /usr/bin/python

import re
import os
import sys
import glob
import subprocess
import csv
from filelock import Timeout, FileLock
import numpy as np


# process private temporary directory
tempDir = "tmp.%d" % os.getpid()

gamma = 0
minK = 4
maxK = 62
model = 'Uniform'

outFile = ''


def loadKmerList( file):
    print("Loading from file %s" % file)
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

    return np.array(seqDict.keys())



def extractKmers( dataset, k, remove = True):
    datasetBasename = os.path.splitext(os.path.basename(dataset))[0]
    kmcOutputPrefix = "k=%d%s" % (k, datasetBasename)
    histFile = "distk=%d_%s.hist" % (k, datasetBasename)
    
    if (not os.path.exists(histFile) or os.path.getsize(histFile) == 0):
        # run kmc on the first sequence
        cmd = "kmc -b -k%d -m2 -fm -ci0 -cs1000000 %s %s %s" % (k, dataset, kmcOutputPrefix, tempDir)
        p = subprocess.Popen(cmd.split())
        p.wait()
        print("cmd: %s returned: %s" % (cmd, p.returncode))
        
        # dump the result -> kmer histogram
        cmd = "kmc_dump %s %s" % ( kmcOutputPrefix, histFile)
        p = subprocess.Popen(cmd.split())
        p.wait()
        print("cmd: %s returned: %s" % (cmd, p.returncode))
    else:
        print("skipping kmer extraction for dataset: %s" % dataset)
        
    # load kmers from histogram file
    vect = loadKmerList(histFile)

    if (remove):
        # remove temporary files
        os.remove(histFile)
        for f in glob.glob(kmcOutputPrefix+'*'):
            os.remove(f)

    return vect




# run jaccard on sequence pair ds with kmer of length = k
def runJaccard( ds1, ds2, k):
    global outFile

    leftKmers = extractKmers(ds1, k, False)
    rightKmers = extractKmers(ds2, k, False)

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
    absentCnt = NMax - (A + B + C) # M01M10M11
    D = absentCnt
    # (M10 + M01) / (M11 + M10 + M01)
    # jaccardIndex = M01M10 / float(M01M10M11)
    jaccardIndex = A / float(NMax - D)
    jaccardDistance = 1 - jaccardIndex
    
    header = ['model', 'gamma', 'k', 'Jaccard Distance', 'A', 'B', 'C', 'D', 'Nmax', 'A/(N-D)']
    data = [model, gamma, k, jaccardDistance, bothCnt, leftCnt, rightCnt, absentCnt, NMax, jaccardIndex]

    lock = FileLock(outFile + '.lck')
    try:
        lock.acquire(5)
        print("Lock acquired.")
        print( data)
        if (not os.path.exists( outFile)):
            f = open(outFile, 'w')
            writer = csv.writer(f)
            writer.writerow(header)
        else:
            f = open(outFile, 'a')
            writer = csv.writer(f)

        writer.writerow(data)
        f.close()
        lock.release()

    except Timeout:
        print("Another instance of this application currently holds the lock.")




def main():
    global minK, maxK, gamma, model, outFile
    
    l = len(sys.argv)
    if (l < 4 or l > 5):
        print("Errore nei parametri:")
        print("Usage: %s seq1 seq2 file.csv [k]" % (sys.argv[0]))
        exit(-1)

    if (l > 3):
        seq1 = sys.argv[1]
        seq2 = sys.argv[2]
        outFile = sys.argv[3]
        
    if (l > 4):
        # solo uno specifico valore di k
        minK = int(sys.argv[4])
        maxK = minK

    # private temporary directory
    tempDir = "tmp.%d" % os.getpid()
    os.mkdir(tempDir)

    m = re.search(r'^(.*)\.*-G=0\.(\d+)', seq2)
    if (m is not None):
        model = m.group(1)
        gamma = float("0."+m.group(2))
    else:
        print('Malformed dataset name')

    print('model: %s, gamma = %.3f' % (model, gamma))

    k = minK
    while (k <= maxK):
        runJaccard( seq1, seq2, k)
        
        if (k < 20):
            k += 2
        elif (k < 30):
            k += 3
        else:
            k += 10

    # cleanup
    os.rmdir(tempDir)



        

if __name__ == "__main__":
    main()
