#! /usr/bin/python

import sys
import os
import os.path
import subprocess


sketchSizes = [ 10000]

class MashData:
    def __init__(self, cmdResults):
        mr = cmdResults.split()
        mashAN = mr[4].decode('UTF-8')
        self.Pv = float(mr[2])
        self.dist = float(mr[3])
        try:
            ns = mashAN.index('/')
            self.A = int(mashAN[:ns])
            self.N = int(mashAN[ns+1:])
        except ValueError:
            self.A = 0
            self.N = 0



def runMash(inputDS1, inputDS2, k):
    # run mash on the same sequence pair
    mashValues = []
    for ss in sketchSizes:
        # extract mash sketch from the first sequence
        cmd = "/usr/local/bin/mash sketch -s %d -k %d %s" % (ss, k, inputDS1)
	print("cmd1: %s" % cmd)      
        p = subprocess.Popen(cmd.split())
        p.wait()

        cmd = "/usr/local/bin/mash sketch -s %d -k %d %s" % (ss, k, inputDS2)
	print("cmd2: %s" % cmd)      
        p = subprocess.Popen(cmd.split())
        p.wait()

        cmd = "/usr/local/bin/mash dist %s.msh %s.msh" % (inputDS1, inputDS2)
	print("cmd3: %s" % cmd)      
        out = subprocess.check_output(cmd.split())

        mashValues.append( MashData( out))

    # dati mash distance
    data2 = []
    for i in range(len(sketchSizes)):
        data2.append( mashValues[i].Pv)
        data2.append( mashValues[i].dist)
        data2.append( mashValues[i].A)
        data2.append( mashValues[i].N)

    # clean up remove kmc temporary files
    os.remove(inputDS1 + '.msh')
    os.remove(inputDS2 + '.msh')

    return data2



def main():
    argNum = len(sys.argv)
    if (argNum != 4):
        print( "Usage:\n%s seq1 seq2 k" % sys.argv[0])
	exit(1)
    else:
        kLen = int(sys.argv[3])
        seq1 = sys.argv[1]
        seq2 = sys.argv[2]

    print(runMash(seq1, seq2, kLen)+"\n")



if __name__ == "__main__":
    main()
