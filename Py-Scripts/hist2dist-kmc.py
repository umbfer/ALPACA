#! /usr/bin/python

import re
import os
import sys
from array import array


inputFile = sys.argv[1]

totalCnt = 0
totalSeqCnt = 0
totalKmer = 0
totalProb = 0.0
totalLines = 0

sumDict = dict()  # dizionario vuoto
dist = 'dist' # directory per le distribuzioni
basename = os.path.splitext(os.path.basename(inputFile))[0]


def main():
    global totalCnt, totalSeqCnt, totalKmer, totalProb, totalLines, sumDict, seqDict

    m = re.search(r'^(.*)_(.*)-(\d+)\.(\d+)', basename)
    if (m is not None):
        model = m.group(2)
        nPairs = int(m.group(3))
        seqLen = int(m.group(4))

    # ogni file contiene l'istogramma di una sola sequenza prodotto con kmc 3
    with open(inputFile) as inFile:
        seqDict = dict()
        totalSeqCnt = 0
        seqNum = 1
        for line in inFile:
            m = re.search(r'^([A-Z]+)[ \t]+(\d+)', line)
            if (m is None):
                print line, " malformed histogram file"
                exit()
            else:
                kmer = m.group(1)
                count = int(m.group(2))

            # totale generale per il calcolo della probabilita' empirica
            totalCnt = totalCnt + count
            totalSeqCnt = totalSeqCnt + count
            totalLines = totalLines + 1

            if kmer in sumDict:
                sumDict[kmer] = sumDict[kmer] + int(count)
            else:
                sumDict[kmer] = int(count)


    print('')

    # freqFile = basename + '.sum'
    # with open(freqFile, "w") as outText:
    #     for key in sumDict:
    #        outText.write("%s\t%d\n" % (key, sumDict[key]))

    totalKmer = 0
    totalProb = 0
    probFile = "%s/%s.dist" % (dist, basename)
    with open(probFile, "w") as outDist :
        kmers = sorted(sumDict.keys())
        # arrLen = [1]
        # arrLen[0] = len(kmers)
        # header = array('i', arrLen )
        # header.tofile(outDist)
        # pv = [0]*arrLen[0]
        # ks = 0
        for key in kmers:
            prob = sumDict[key] / float(totalCnt)
            # pv[ks] = prob
            # ks = ks + 1
            totalProb = totalProb + prob
            totalKmer = totalKmer + 1
            outDist.write("%s\t%.10f\n" % (key, prob))

        # salva l'intero vettore di probabilita'
        # fa = array('d', pv)
        # fa.tofile(outDist)

    print("total kmer values:\t%d" % totalLines)  # numero dei conteggi
    print("total distinct kmers:\t%d" % totalKmer)  # numero kmers
    print("total kmers counter:\t%d" % totalCnt)  # totale conteggio
    print("total prob-distr.:\t%f" % totalProb)  # totale distribuzione di probabilita'


if __name__ == "__main__":
    main()
