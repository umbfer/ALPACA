#! /usr/bin/python

import re
import os
import sys

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


# salva la distribuzione ed istogramma della sequenza prevSeq
def SaveSeqFiles(prevSeq, model, len):
    pairSeq = 'A' if (prevSeq % 2 == 0) else 'B'
    d1 = "%s/%s" % (seqDistDir, model)
    d2 = "%s/%d" % (d1, len)
    if not os.path.exists(d1):
        os.mkdir(d1)
    if not os.path.exists(d2):
        os.mkdir(d2)
    fileName1 = "%s/%s-%s-%s.dist" % (d2, basename, str((prevSeq / 2) + 1).zfill(5), pairSeq)
    fileName2 = "%s/%s-%s-%s.hist" % (d2, basename, str((prevSeq / 2) + 1).zfill(5), pairSeq)
    with open(fileName1, "w") as outFileDist:
    #    with open(fileName2, "w") as outFileHist:
        for key in sorted(seqDict.keys()):
            prob = seqDict[key] / float(totalSeqCnt)
            outFileDist.write("%s\t%.10f\n" % (key, prob))
            # poi salva l'istogramma della sequenza precedente
            # outFileHist.write("%s\t%d\n" % (key, seqDict[key]))


def main():
    global totalCnt, totalSeqCnt, totalKmer, totalProb, totalLines, writeSequenceHistogram, sumDict, seqDict

    prevSeq = -1

    if (writeSequenceHistogram):
        if not os.path.exists(seqDistDir):
            os.mkdir(seqDistDir)

    m = re.search(r'^(.*)_(.*)-(\d+)\.(\d+)', basename)
    if (m is None):
        print basename, " malformed histogram filename"
        exit()
    else:
        model = m.group(2)
        nPairs = int(m.group(3))
        len = int(m.group(4))

    with open(inputFile) as inFile:
        for line in inFile:
            m = re.search(r'^\((\d+),\(([A-Z]+),(\d+)\)\)', line)
            if (m is None):
                print line, " malformed histogram file"
                exit()
            else:
                seqNum = int(m.group(1))
                kmer = m.group(2)
                count = int(m.group(3))

            if (seqNum != prevSeq):
                if (writeSequenceHistogram and prevSeq != -1):
                    # salva la distribuzione e l'istogramma della sequenza precedente
                    SaveSeqFiles(prevSeq, model, len)

                # abbiamo salvato prevSeq prossima sequenza seqNum
                seqDict = dict()
                totalSeqCnt = 0
                prevSeq = seqNum
                sys.stdout.write('\r%d / %d Complete\r' % (seqNum,nPairs*2)),
                sys.stdout.flush()

            # totale generale per il calcolo della probabilita' empirica
            totalCnt = totalCnt + count
            totalSeqCnt = totalSeqCnt + count
            totalLines = totalLines + 1

            if kmer in sumDict:
                sumDict[kmer] = sumDict[kmer] + int(count)
            else:
                sumDict[kmer] = int(count)

            if kmer in seqDict:
                seqDict[kmer] = seqDict[kmer] + int(count)
            else:
                seqDict[kmer] = int(count)

    if (writeSequenceHistogram and prevSeq != -1):
        # salva la distribuzione e l'istogramma dell'ultima sequenza
        SaveSeqFiles(seqNum, model, len)

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
