#! /usr/local/bin/python

import re
import os
import sys

writeSequenceHistogram = True
seqDistDir = 'seqDists'



def main():

    CheckDistribution(sys.argv[1])

    
    
def CheckDistribution(inputFile):

    seq = 0
    (Acnts,Ccnts,Gcnts,Tcnts,Wcnts,Tot) = [0,0,0,0,0,0]
    seqName = ''
    with open(inputFile) as inFile:
        for line in inFile:
            if (line.startswith(">")):
                seq += 1
                if (seq > 1):
                    print("seq %s:\tA=%d (%.2f%%), C=%d (%.2f%%), G=%d (%.2f%%), T=%d (%.2f%%), bad=%d, tot=%d\n" %
                          (seqName, Acnts, Acnts/Tot, Ccnts, Ccnts/Tot, Gcnts, Gcnts/Tot, Tcnts, Tcnts/Tot, Wcnts, Tot))

                seqName = line[1:]
                (Acnts,Ccnts,Gcnts,Tcnts,Wcnts,Tot) = [0,0,0,0,0,0]
            else:
                for b in line:
                    if (b == '\n'):
                        continue
                    b = b.upper()
                    Tot += 1
                    match b:
                        case 'A':
                            Acnts += 1
                        case 'C':
                            Ccnts += 1
                        case 'G':
                            Gcnts += 1
                        case 'T':
                            Tcnts += 1
                        case _:
                            Wcnts += 1
                            # print("bad = %d\n" % ord(b))

        print("seq %s:\tA=%d (%.2f%%), C=%d (%.2f%%), G=%d (%.2f%%), T=%d (%.2f%%), bad=%d, tot=%d\n" %
                          (seqName, Acnts, Acnts/Tot, Ccnts, Ccnts/Tot, Gcnts, Gcnts/Tot, Tcnts, Tcnts/Tot, Wcnts, Tot))




if __name__ == "__main__":
    main()
