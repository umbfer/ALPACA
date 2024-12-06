#! /usr/bin/python

import re
import os
import sys
import time
from array import array

inputFile = sys.argv[1]

totalCnt = 0
totalSeqCnt = 0
totalKmer = 0
totalProb = 0.0
totalLines = 0
writeSequenceHistogram = True

sumDict = dict()  # dizionario vuoto
seqDict = dict()

seqDistDir = 'seqDists'
basename = os.path.splitext(os.path.basename(inputFile))[0]


def OutputFileName():
    return( "%s/%s-All.distbin" % (seqDistDir, basename))


# salva la distribuzione ed istogramma della sequenza prevSeq
def SaveDistributions(model, pairs):
    fileName1 = OutputFileName()
    with open(fileName1, "w") as outFileDist:
        kmers = sorted(seqDict.keys())
        arrLen = [1]
        arrLen[0] = len(kmers)
        header = array('i', arrLen )
        header.tofile(outFileDist)
        pv = [0]*arrLen[0]
        for seq in range(2 * pairs):
            ks = 0
            for key in sorted(seqDict.keys()):
                av = seqDict[key]
                # outFileDist.write("%s " % (key))
                pv[ks] = av[seq] / float(totalSeqCnt)
                ks = ks + 1

            # salva l'intero vettore di probabilita'
            fa = array('d', pv)
            fa.tofile(outFileDist)
            #outFileDist.write("\n")



def main():
    global totalCnt, totalSeqCnt, totalKmer, totalProb, totalLines, writeSequenceHistogram, sumDict, seqDict

    prevSeq = 0 # N.B. la  prima sequenza ha indice 0

    if (writeSequenceHistogram):
        if not os.path.exists(seqDistDir):
            os.mkdir(seqDistDir)

    m = re.search(r'^(.*)k=(\d+)_(.*)-(\d+)\.(\d+)', basename)
    if (m is None):
        print basename, " malformed histogram filename"
        exit()
    else:
        kLen = int(m.group(2))
        model = m.group(3)
        nPairs = int(m.group(4))
        len = int(m.group(5))

    # crea il file distribuzione solo se l'input e' piu' recente
    mt1 = os.path.getmtime(inputFile)
    out = OutputFileName()
    if (os.path.exists(out)):
        mt2 = os.path.getmtime(out)
        if (mt2 > mt1):
            print( "The distribution file already exists and it is more recent than the histogram file")
            print( "%s -> %s vs\n%s -> %s" % (inputFile, time.ctime(mt1), out, time.ctime(mt2)))
            exit(-1)

    with open(inputFile) as inFile:
        for line in inFile:
            # pattern ( idSeq, ( k-mer, cnt))
            m = re.search(r'^\((\d+),\(([A-Z]+),(\d+)\)\)', line)
            if (m is None):
                print line, " malformed histogram file"
                exit(-1)
            else:
                seqNum = int(m.group(1))
                kmer = m.group(2)
                count = int(m.group(3))

            if (seqNum != prevSeq):
                # salviamo tutto il dizionario unico alla fine
                if (totalSeqCnt != len - kLen + 1):
                    print("Wrong number of k-mers %d vs %d" & (totalSeqCnt, len - kLen + 1))
                    exit(-1)
                totalSeqCnt = 0
                prevSeq = seqNum
                sys.stdout.write('\r%d / %d Complete\r' % (seqNum,nPairs*2)),
                sys.stdout.flush()

            # totale generale per il calcolo della probabilita' empirica
            totalCnt = totalCnt + count
            totalSeqCnt = totalSeqCnt + count
            totalLines = totalLines + 1

            if kmer in sumDict:
                sumDict[kmer] = sumDict[kmer] + count
            else:
                sumDict[kmer] = count

            if kmer in seqDict:
                av = seqDict[kmer]
                if (av[seqNum] != 0):
                    print("Errore seqNum %d non zero initial value" & seqNum)

                av[seqNum] = count #otherwise av[seqNum] + count ????
                seqDict[kmer] = av
            else:
                av = [0]*(2*nPairs)
                av[seqNum] = count
                seqDict[kmer] = av

    if (writeSequenceHistogram):
        # salva la distribuzione e l'istogramma dell'ultima sequenza
        SaveDistributions(model, nPairs)

    print('')

    freqFile = basename + '.sum'
    with open(freqFile, "w") as outText:
        for key in sumDict:
            outText.write("%s\t%d\n" % (key, sumDict[key]))

    totalKmer = 0
    totalProb = 0
    probFile = basename + '.dist'
    with open(probFile, "w") as outText:
        for key in sumDict:
            prob = sumDict[key] / float(totalCnt)
            totalProb = totalProb + prob
            totalKmer = totalKmer + 1
            outText.write("%f\n" % (prob))

    print("total sequence number:\t%d" % (int(seqNum) + 1))  # numero sequenze seqId starts from 0
    print("total kmer values:\t%d" % totalLines)  # numero dei conteggi
    print("total distinct kmers:\t%d" % totalKmer)  # numero kmers
    print("total kmers counter:\t%d" % totalCnt)  # totale conteggio
    print("total prob-distr.:\t%f" % totalProb)  # totale distribuzione di probabilita'


if __name__ == "__main__":
    main()
