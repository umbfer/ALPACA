#! /usr/local/bin/python

import re
import os
import sys
import random




def hamming_distance(seq1: str, seq2: str) -> int:
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    d = sum([ c1 != c2 for c1, c2 in zip(seq1, seq2)])
    return d

def hamming_distance2(seq1: str, seq2: str) -> int:
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(seq1, seq2))))



# parametri sulla linea di comando
# inputSeqence theta
def CompareSequences():

    (dist, tot1, tot2) = (0, 0, 0)

    if (len(sys.argv) == 3):
        inputFile1 = sys.argv[1]
        inputFile2 = sys.argv[2]
    else:
        print("Errore nei parametri:\nUsage: %s sequence1 sequence2" % os.path.basename(sys.argv[0]))
        exit(-1)

    with open(inputFile1, "r") as inFile1, open(inputFile2, "r") as inFile2:
        while( True):
            # skip all comment lines from file1
            while (True):
                line1 = inFile1.readline()
                if (not line1.startswith(">")):
                    break

            while (True):
                line2 = inFile2.readline()
                if (not line2.startswith(">")):
                    break
            if (line1 == "" and line2 == ""):
                # entrambi i file sono terminati
                break

            dist += hamming_distance( line1, line2)
            tot1 += len(line1) - 1 # \n
            tot2 += len(line2) - 1 # \n

        inFile2.close()
    print("Hamming distance: %s (%d) vs %s (%d) = %d (%.2f%%)" % (inputFile1, tot1, inputFile2, tot2, dist, dist/tot1*100))



if __name__ == "__main__":
    CompareSequences()

