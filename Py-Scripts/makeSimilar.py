#! /usr/local/bin/python

import re
import os
import sys

writeSequenceHistogram = True
seqDistDir = 'seqDists'

# parametri sulla linea di comando
# seqenceName len1 base1 len2 base2
def ModifySequence():

    basis = []
    if (len(sys.argv) == 5):
        inputFile = sys.argv[1]
        newLen = int(sys.argv[2])
        basis.append( sys.argv[3])
        basis.append( sys.argv[4])
    else:
        print("Errore nei parametri:\nUsage: %s sequenceName #lenToAdd base1 base2 ")
        exit(-1)

    outFile = "%s-mod.fasta" % (os.path.splitext( inputFile)[0])
    seq = 0
    newSeq = []
    with open(inputFile) as inFile:
        for line in inFile:
            if (line.startswith(">")):
                newSeq.append(line)
            else:
                seq1 = line[0:-1] + basis[seq] * (newLen / len(basis[seq]))
                newSeq.append(seq1)
                seq += 1

    with open(outFile, "w") as outText:
        outText.write("%s%s\n" % (newSeq[0], newSeq[1])) # \n are in the original strings
        outText.write("%s%s\n" % (newSeq[2], newSeq[3])) # \n are in the original strings

    print("%s -> %s" % (inputFile, outFile))

if __name__ == "__main__":
    ModifySequence()
