#! /usr/bin/python

import re
import os
import sys
import shutil

rootDir = sys.argv[1] if (len(sys.argv) > 1) else '/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/dataset5-1000'

mashDir = 'mash'
mergeDir = 'merge'
kValues = [ 4, 6, 8, 10]

mash = dict()  # dizionario vuoto
other = dict()


def main():
    global mashDir, mergeDir, rootDir, kValues

    fileCnt = 0
    for k in kValues:
        # cancella la directory destinazione
        mergePath = '%s/%s/k=%d' % ( rootDir, mergeDir, k)
        if os.path.exists(mergePath):
            shutil.rmtree(mergePath)

        os.mkdir(mergePath)

        mashPath = "%s/%s/k=%d" %( rootDir, mashDir, k)
        for f in os.listdir( mashPath):
            if (f.endswith(".csv")):
                file = '%s/%s' % (mashPath, f)
                print("Processing: %s" % file)
                fileCnt = fileCnt + 1
            else:
                continue

            mash = dict()
            cnt = 0
            # carica il file risultati (mash) da fondere
            with open(file) as inFile:
                for line in inFile:
                    cnt = cnt + 1
                    ll = line.split(',')
                    if (cnt == 1):
                        measureName = ll[2].rstrip()
                    else:
                        key = "%s,%s" % (ll[0], ll[1])
                        mash[ key] = ll[2].rstrip()

            if (cnt != 1001):
                print("Only %d results processed from file %s" % (cnt, file))
            
            otherPath = "%s/k=%d/%s" %( rootDir, k, f)
            other = dict()
            cnt = 0
            with open(otherPath) as inFile:
                for line in inFile:
                    cnt = cnt + 1
                    ll = line.split(',')
                    if (cnt == 1):
                        measuresHdr = line.rstrip()
                    else:
                        key = "%s,%s" % (ll[0], ll[1])
                        other[ key] = line.rstrip() # riuove il newline alla fine

            if (cnt != 1001):
                print("Only %d results processed from file %s" % (cnt, otherPath))

            # salva il file merged
            mergePath = '%s/%s/k=%d/%s' % ( rootDir, mergeDir, k, f)
            with open(mergePath, "w") as outFileCSV:
                outFileCSV.write("%s,%s\n" % (measuresHdr, measureName))
                for key in sorted(other.keys()):
                    if (mash.has_key(key)):
                        outFileCSV.write("%s,%s\n" % (other[key], mash[key]))
                    else:
                        print("chiave %s: mancante nel file %s" % (key, f))
                        sys.exit("file mash incompleto")

    print("Processed %d file" % fileCnt);



if __name__ == "__main__":
    main()
