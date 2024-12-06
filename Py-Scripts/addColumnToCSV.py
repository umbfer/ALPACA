#! /usr/local/bin/python3

import re
import os
import sys
import csv


def main():
    if (len(sys.argv) < 2):
        print("Wrong number of arguments. Usage:\n%s k-value file1 [file2 ...]", sys.argv[0])
        exit(-1)

    for file in sys.argv[2:]:
        addKColumn(int(sys.argv[1]), file)

    
    
def addKColumn( k, inputFile):

    # dist-k=4_MotifRepl-Sh-1000.10000.G=010.csv or dist-k=4_ShuffledEColi-1000.10000.csv
    m = re.search(r'^dist-k=(\d+)_(.*)\.csv', inputFile)
    if (m is None):
        print("Malformed file name %s" % inputFile)
        exit(-1)
    else:
        if (k != int(m.group(1))):
            print("Wrong k value: %d vs %s" % (k, inputFile))
            exit(-2)

        outfileName = m.group(2) + '.csv'
        outfileName = '../fadeResults-k=%d.csv' % k
        m = re.search(r'^dist-k=(\d+)_(.*)-1000\.(\d+)(.*)\.csv', inputFile)
        if (m is None):
            print("Malformed file name %s" % inputFile)
            exit(-1)
        else:
            model = m.group(2)
            seqLen = int(m.group(3))
            gamma = 0. if (len( m.group(4)) == 0) else float("0." + m.group(4)[3:])
            print( "model = %s, seqLen = %d, k = %d, gamma = %f" % (model, seqLen, k, gamma))

    header = ['Model', 'gamma', 'seqLen', 'seqId', 'k', 'D2', 'D2z', 'Euclidean']
    l = 1
    with open(outfileName, 'a', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)

        with open(inputFile) as inFile:
            for line in inFile:
                if (l > 1):
                    line = line.strip()
                    vals = line.split(',')
                    # MotifRepl.00087.G=010-A
                    m = re.search(r'^(.+)\.(\d+)(.*)', vals[0])
                    if (m is None):
                        print("regular expression failed: %s" % vals[0])
                        exit(-3)
                    else:
                        # if (model != m.group(1)):
                        #    print("wrong model %s vs %s" % ( model, m.group(1)))
                        #    exit(-4)
                        # else:
                        seqId = int(m.group(2))

                    # write a row
                    data = [ model, gamma, seqLen, seqId, k, vals[2], vals[3], vals[4]]
                    writer.writerow(data)
                else:
                    l = l + 1



if __name__ == "__main__":
    main()
