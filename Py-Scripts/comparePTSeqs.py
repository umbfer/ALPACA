#! /usr/bin/python

import re
import os
import sys
import subprocess


inputFile1 = sys.argv[1]
inputFile2 = sys.argv[2]

totalCnt = 0
totalSeqCnt = 0
totalKmer = 0
totalProb = 0.0
totalLines = 0
writeSequenceHistogram = True

sumDict = dict()  # dizionario vuoto
seqDict = dict()

seqDistDir = 'seqDists'
#basename = os.path.splitext(os.path.basename(inputFile))[0]




def main():
    global totalCnt, totalSeqCnt, totalKmer, totalProb, totalLines, writeSequenceHistogram, sumDict, seqDict

    p = subprocess.Popen(['cmp', '-l', inputFile1, inputFile2],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    skip = 0
    offset = skip
    cnt = 0
    maxLen = 0

    for line in out.splitlines():
        m = re.search(r'(\d+)[ \t]+(\d+)[ \t]+(\d+)', line)
        if (m is None):
            print line, " unknown input format"
            exit()
        else:
            actualOffset = int(m.group(1))

        if ((actualOffset - offset) >= 10):
            cnt = cnt + 1
            len = actualOffset - offset
            if (len > maxLen):
                maxLen = len
            print("found Pattern Transferred @ %d, len = %d, cnt = %d" % (offset - skip, len, cnt))
            if (cnt > 10):
                break

        offset = actualOffset

    print("%d transferts, maxLength = %d" % (cnt, maxLen))


if __name__ == "__main__":
    main()
