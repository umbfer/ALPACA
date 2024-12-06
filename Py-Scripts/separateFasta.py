#! /usr/bin/python

import re
import os
import sys
import random
import numpy as np

#
# Usage:
# separateFasta.py sequence.fasta gamma
#
# Crea una nuova sequenza sequence-G=gamma.fasta in cui ogni base con probabilita' gamma
# viene rimpiazzata (da una delle altre 3)
#
if (len(sys.argv) != 3):
    print("Errore nei parametri.Usage:\n%s inputSequence.fasta gamma" % sys.argv[0])
    exit(-1)


inputFile = sys.argv[1]
gamma = float(sys.argv[2])
outFile = '%s-G=%.3f.fasta' % (os.path.splitext(os.path.basename(inputFile))[0], gamma)
probRange = 10000

bases = ['A', 'C', 'G', 'T']

# seleziona a caso uno degli altri 3 caratteri
def randomBase( base):
    b = base.upper()
    l = []
    for c in bases:
        if ( b != c):
            l.append(c)

    return l[random.randint(0, 2)]


def main():

    cnt = 0
    totCnt = 0
    out = open(outFile, "w")
    with open(inputFile) as inFile:
        for line in inFile:
            if (line.startswith('>')):
                out.write( line.rstrip('\n') + ' G=%.3f\n' % gamma)
            else:
                ll = len(line)
                totCnt += ll
                nl = np.empty( ll, np.character)
                for i in range( ll):
                    rn = random.randint(0, probRange)
                    c = line[i]
                    if (c in bases and rn < gamma * probRange):
                        # ok substitute c with another base randomly choosed
                        nl[i] = randomBase( c)
                        cnt += 1
                    else:
                        # otherwise keep the previous char
                        nl[i] = c

                out.write( nl) # \n are in the original strings

                # sys.stdout.write('.')
                # sys.stdout.flush()

    out.close()
    print('terminated, %d / %d bases substituted' % (cnt, totCnt))



if __name__ == "__main__":
    main()
