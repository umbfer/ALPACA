#! /usr/local/bin/python3

import re
import os
import sys
import string
import random
import time
import csv
from pathlib import Path
import numpy as np


ext = '.fna'
basis = string.ascii_uppercase + string.digits + ' '

# parametri sulla linea di comando
# inputSeqence theta

def main():
    k = 0
    if (len(sys.argv) == 3):
        inputPath = sys.argv[1]
        k = int(sys.argv[2])
    else:
        print("Errore nei parametri:\nUsage: %s InputSequence kmerLength" % os.path.basename(sys.argv[0]))
        exit(-1)

    inputFile = os.path.basename(inputPath)
    baseName, ext = os.path.splitext( inputPath)
    outFile = f"{baseName}-{int(time.time())}.csv"
    f = open(outFile, 'w')
    writer = csv.writer(f)
    header = ['Input File', 'k', 'theta', 'A', 'B', 'C', 'D', 'N', 'Jaccard',
              'Total Distinct left', 'Total Kmers left', 'Total Distinct right', 'Total Kmers right']
    writer.writerow(header)

    for k in range(2,11,1):
        (totalLeft, hist1) = kmerExtraction(inputPath, k, 0)
        leftKeys =  np.array(list(hist1.keys()))

        for theta in [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]:
            # inputFile2 = MoveAwaySequence(inputPath, theta)
            (totalRight, hist2) = kmerExtraction(inputPath, k, theta)

            # leftKeys =  np.array(list(hist1.keys()))
            rightKeys =  np.array(list(hist2.keys()))

            intersection = np.intersect1d( leftKeys, rightKeys)
            A = bothCnt = intersection.size
            B = leftCnt = leftKeys.size - bothCnt
            C = rightCnt = rightKeys.size - bothCnt

            NMax = pow(len(basis), k)
            D = absentCnt = NMax - (A + B + C) # NMax - M01M10M11

            # Jaccard dissimilarity => Jaccard = 1 - A/(N - D)
            try:
                jaccard = 1 - A / (NMax - D)
            except ZeroDivisionError:
                jaccard = 'UnDefined'

            print(f"Distance: {jaccard} k: {k} theta: {theta}, A: {A:,}, B: {B:,}, C: {C:,}, D: {D:,}")

            writer.writerow( [ inputFile, k, theta, A, B, C, D, NMax, jaccard,
                              leftKeys.size, totalLeft, rightKeys.size, totalRight])

    f.close()




def kmerExtraction(inputPath: string, k: int, theta: int):

    inputFile = Path(inputPath).name
    print( "*********************************************************")
    print( "Extracting kmers from sequence: %s for k: %d and theta: %d" % (inputFile, k, theta))

    count = {}
    (cnt, subst, totalLen) = (0, 0, 0)
    regex0 = re.compile( '[\r\n]')
    regex1 = re.compile( '[%s]' % re.escape(string.punctuation))
    regex2 = re.compile('\s\s+')
    regex3 = re.compile('[\u0080-\uFFFF]')

    allText = ""
    with open(inputPath, encoding='latin-1') as inFile:
        for line in inFile:
            # 1) trasforma tutti i caratteri alfanumeri in uppercase
            l0 = line.upper()
            # 2) sostituisce newline con spazio
            l1 = regex0.sub(' ', l0)
            # 3) rimpiazza tutti i caratteri di punteggiatura con spazio
            l2 = regex1.sub(' ', l1)
            # 4) rimpiazza più spazi consecutivi in un unico spazio
            l3 = regex2.sub(' ', l2)
            # 5) rimuove i caratteri ascii > 127
            line = regex3.sub( '', l3)

            if theta > 0:
                s = list(line)
                for i in range(len(line)):
                    if (random.randrange(100) < theta):
                        # l'elemento i-esimo viene sostituito
                        b = s[i]
                        newBase = b
                        while( newBase == b):
                            newBase = random.choice( basis)
                        s[i] = newBase
                        subst += 1
                line = "".join(s)

            # e concatena TUTTO l'input in una unica stringa
            allText = allText + line
            cnt += 1
            # print(i, l3)
            totalLen += len(line)

    print(f"{inputFile} -> {subst:,}/{totalLen:,} substitutions")

    totalLen = len(allText)
    last = totalLen - k + 1

    for i in range(last):
        kmer = allText[i:i+k]
        c = count.get(kmer, 0)
        count[kmer] = c+1

    # [print(f"{key}: {value}") for key, value in count.items() if (value > 1)]

    print(f"Total Distinct: {len(count.keys()):,} / Total: {last:,} kmers (over {len(basis)**k:,} possible kmers)")
    return (totalLen, count)



def MoveAwaySequence(inputFile: string, theta: int):

    baseName, ext = os.path.splitext( inputFile)
    outFile = "%s-%02d%s" % (baseName, theta, ext)
    if (os.path.exists(outFile)):
        os.remove(outFile)
        # print("Output File: %s already exists. Exiting." % outFile)
        # return outFile

    print( "*********************************************************")
    print( "Creating sequence: %s from sequence: %s theta: %d" % (Path(outFile).name, Path(inputFile).name, theta))
    print( "*********************************************************")

    (written, subst, totLen) = (0, 0, 0)
    newBase = ''
    out = []

    with open(outFile, "w") as outText:
        with open(inputFile, encoding='latin-1') as inFile:
            for line in inFile:
                s = list(line)
                for i in range(len(line)):
                    if (random.randrange(100) < theta):
                        # l'elemento i-esimo viene sostituito
                        b = s[i]
                        newBase = b
                        while( newBase == b):
                            newBase = random.choice( basis)
                        s[i] = newBase
                        subst += 1

                out = "".join(s)
                totLen += len(line) - 1

                outText.write(out) # \n are in the original strings
                m = totLen // 1048576

    print(f"{outFile} -> {subst:,}/{totLen:,} substitutions")
    return outFile





if __name__ == "__main__":
    main()
