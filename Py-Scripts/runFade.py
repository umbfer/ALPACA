#!/usr/bin/python

#
# $Id: runFade.py 1738 2020-12-11 12:25:28Z cattaneo@dia.unisa.it $
#

import sys
import datetime
import re
import math
import ToolBox
from os.path import basename, splitext
from optparse import OptionParser

Version = "$Revision: 1738 $".split(' ')[1]


# numero executors
deployMode = 'client' # client / cluster (immagino che client mode liberi uno slave, facendo girare l'AM sul master)
numExec = 48
numBins = 267
kVal = 0
gVal = 0
masterConf = "confs/allMeasures.conf"
confFile   = "fade.conf"
logFile    = "run.log"
inputPath  = "data/dataset5-1000"
resultsPrefix = "results"
fadePath = " build/libs"

kValues = (4, 6, 8, 10)

def ComputeKValue( seqLen):
    # somma(len(i)) / n
    return math.ceil(math.log( float(seqLen), 4)) - 1
    


def is_int(val):
    try:
        num = int(val)
    except ValueError:
        return False
    return True
                

class Logger(object):

    def __init__(self, scriptName):
        self.terminal = sys.stdout
        ct = datetime.datetime.now()
        filename = "%s-%s.log" % (scriptName, ct.strftime("%Y%m%d-%H%M%S"))
        self.log = open(filename, "a")
        self.fileno = sys.stdout.fileno

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        pass    

    def fileno(self):
        #this flush method is needed for subprocess  compatibility.
        return( self.fileno)


    

def main():

    global options, kValues, numExec
        
    ScriptName = basename(splitext(__file__)[0])

    sys.stdout = Logger( ScriptName)
    
    parser = OptionParser(usage="usage: %prog [options]",
                          version='%s %s' % (ScriptName, Version))
    
    parser.add_option("-x", "--executors", action="store",
                      dest="executorsNumber", default = str(numExec),
                      help="Number of spark task executors")
    
    
    parser.add_option("-f", "--nokspec", action="store_true", # optional
                      dest = "nokspec", default = False,
                      help="Avoid the script to compute and specify the k-value")
    
    parser.add_option("-k", "--kvalues", action="append", # optional
                      dest = "clkvalues", default = [],
                      help="Specify one or more alternate values for k to replace the one automatically computed depending on the sequence length")
    
    parser.add_option("-b", "--bins", action="store", # optional
                      dest = "binNumber", default = "",
                      help="Overide the default number of bins used for partial aggregation")
    
    parser.add_option("-d", "--debug", action="store_true",
                      dest="debug_flag", default = False,
                      help="Print debug information and do not execute commands")
    
    parser.add_option("-i", "--inputDataSet", action="store",
                      dest="inputDataSet", default = inputPath,
                      help="Set the full HDFS path of directory contaoining the input dataset(s). It may be a regular expression to select only some ds")

    parser.add_option("-p", "--fadePath", action="store",
                      dest="fade", default = fadePath,
                      help="Set the path of the directory containing the FADE fade-1.0.0-all.jar file.")
    
    parser.add_option("-l", "--logFile", action="store", # ==> optional
                      dest="logFile", default=logFile,
                      help="Specify log filename")

    parser.add_option("-c", "--confFile", action="store", # ==> optional
                      dest="confFile", default=masterConf,
                      help="Specify master configuration  filename")
    
    (options, args) = parser.parse_args()

    numExec = int(options.executorsNumber)
    if (len(options.clkvalues) > 0):
        kValues = map(int, options.clkvalues) 
    
    dir = ToolBox.GetCommandOutput("hdfs dfs -ls %s" % options.inputDataSet)

    if (not dir):
        print( 'No input file specifiend in the directory or pattern %s' % (options.inputDataSet))
        exit(0)
    else:
        dirOut = dir.split('\n')

    s = 1 if (dirOut[0].startswith('Found')) else 0
    dsList = []
    
    for f in dirOut[s:]: 
    	dsList.append(f.split()[7])
	

    ToolBox.RunCommand( "rm %s" % options.logFile)

    resultsPath = "%s/%s" % (resultsPrefix, basename(splitext(options.confFile)[0]))

    ToolBox.RunCommand( "hdfs dfs -mkdir %s" % (resultsPath))

    print( "Test starting, script %s version: %s, %d dataset selected" % (ScriptName, Version, len(dsList)))
           
    for ds in dsList:

        m = re.search( r'.*-(\d+)\.(\d+)\..*', ds)
        if (m is None):
            print(ds+ " malformed input dataset")
            exit()
        else:
            seqNum = int(m.group(1))
            seqLen = int(m.group(2))

        m = re.search( r'.*\.G=0\.(\d+)\.fasta', ds)
        gVal = 0 if (m is None) else int(m.group(1))

        for kVal in kValues:

            # kVal   = 0 if (options.nokspec) else ComputeKValue( seqLen) if (options.clkvalues == "") else int(options.clkvalues)
 
            if (len(options.binNumber) > 0 and is_int(options.binNumber) and
                int(options.binNumber) > 0):
                slices = int(options.binNumber)
            else:
                if (kVal <= 6):
                    slices = 73
                elif (kVal <= 8):
                    slices = 269
                elif (kVal == 9):
                    slices = 1023
                elif (kVal == 10):
                    slices = 2069
                elif (kVal == 11):
                    slices = 4079
                else:
                    slices = 8069

            inputDataSet = splitext(basename( ds))[0]

            P1 = "k=" if (options.nokspec) else "k=%d" % (kVal) # per spaced words is null
            P2 = "output=%s/k=%d/dist-%s_%s" % (resultsPath, kVal, P1, inputDataSet)
            P3 = "slices=%d" % slices
            P4 = "b=%d" % slices

            datasetType = inputDataSet[0:inputDataSet.rindex('-')] # at least one dash must be there

            with open( options.confFile, "r") as sources:
                lines = sources.readlines()

            with open( confFile, "w") as config:
                for line in lines:
                    config.write(re.sub( r"^input=", r"input=" + ds,
                                 re.sub( r"^k=", P1,
                                         re.sub( r"^b=", P4,
                                                 re.sub( r"^slices=", P3,
                                                             re.sub( r"^output=", P2, line))))))

            # Eseguo fade
            dtStart = datetime.datetime.now()
            out = open( logFile, 'a')

            runSpark = "spark-submit --master yarn --deploy-mode %s \
                --driver-memory 22g --num-executors %d \
                --executor-cores 7 --executor-memory 27500m \
                --conf spark.eventLog.enabled=true\
                --conf spark.hadoop.dfs.replication=1\
                --conf spark.yarn.executor.memoryOverhead=2500m\
                %s/fade-1.0.0-all.jar %s" % (deployMode, numExec, options.fade, confFile)

            print( "%s mode: %s Dataset: %s, pairs: %d, len: %d, k: %d, b: %d, g: %d" %
                   (dtStart.strftime("%Y/%m/%d, %H:%M:%S"), deployMode,
                    datasetType, seqNum, seqLen, kVal, slices, gVal))
            if (not options.debug_flag):
                retCode = ToolBox.RunCommand2(runSpark, out)
            else:
                print(runSpark)

            out.close()
            # print  runSpark

            dtEnd = datetime.datetime.now()
            elapsed = (dtEnd - dtStart).seconds

        	# raw_input("Press the <ENTER> key to continue...")

	 



if __name__ == "__main__":
        main()
