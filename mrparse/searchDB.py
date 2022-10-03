#! /usr/bin/env ccp4-python
#
#     Copyright (C) 2022 Ronan Keegan
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Application.
#     A copy of the CCP4 licence can be obtained by writing to the
#     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#
# A wrapper for phmmer
#

import os, sys
import subprocess
import shlex
import time
import pickle

from mrbump.seq_align import simpleSeqID
from mrbump.tools import MRBUMP_utils


class PHHit:
    def __init__(self):
        self.chainName = ""
        self.afdbName = ""
        self.chainID = ""
        self.domainID = ""
        self.rank = ""
        self.prob = 0.0
        self.evalue = 0.0
        self.pvalue = 0.0
        self.score = 0.0
        self.ndomains = 0
        self.domScores = dict([])
        self.alignment = ""
        #self.alnRange = ""
        self.alnRange = [0,0]
        self.ecodRange = []
        self.tarRange = [0,0]
        self.tarExtent = 0
        self.tarMidpoint = 0.0
        self.cols = 0
        self.localSEQID = 0.0
        self.overallSEQID = 0.0
        self.resolution = 0.0
        self.expdta = None
        self.releaseDate = None
        self.alignments = dict([])
        self.modelResRange = None
        self.modelResStart = None
        self.modelResEnd   = None

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


class Domains:
    def __init__(self):
        self.targetName = ""
        self.ID = 0
        self.midpoint = 0
        self.extent = 0
        self.matches = []
        self.ranges = []

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


class phmmer:
    def __init__(self):
        self.phmmerEXE = os.path.join(os.environ["CCP4"], "libexec", "phmmer")
        self.cpu = 1
        self.seqin = ""
        self.seqdb = ""
        self.seqout = ""
        self.logfile = ""
        self.alnfile = ""
        self.iterations = 1
        self.termination = False
        self.workingDIR = ""

        self.extentTolerance = 50
        self.midpointTolerance = 20

        self.resultsList = []
        self.resultsDict = dict([])
        self.targetDomainDict = dict([])

    def getAlignments(self):
        """ Get the alignments for each of the matches in the results list """

        if os.path.isfile(self.alnfile) == False:
            sys.stdout.write("Phmmer Error: can't find alignments file: \n  %s\n" % self.alnfile)
            sys.exit()

        f = open(self.alnfile, "r")
        lines = f.readlines()
        f.close()

        for line in lines:
            for i in self.resultsList:
                if i[0:6] in line and "#=" not in line[0:2] and (line.split()[0].split("/")[-1] == self.resultsDict[i].alnRange):
                    self.resultsDict[i].alignment += line.split()[-1]
                    self.resultsDict[i].alnRange = line.split("/")[-1].split()[0]

    def getPhmmerAlignments(self, targetSequence="", phmmerALNLog="", PDBLOCAL=None, DB=None, seqMetaDB=None):
        """ Extract the alignments from the Phmmer logfile """

        if os.path.isfile(self.logfile) == False:
            sys.stdout.write("Phmmer Error: can't find log file: \n  %s\n" % self.logfile)
            sys.exit()

        # Collect other details from the run log
        CAPTURE = False
        count = 1
        scoreList=[]
        TEMPresultsDict = dict([])
        for line in phmmerALNLog:

            if "No hits satisfy inclusion thresholds; no alignment saved" in line:
                sys.stdout.write("Sorry, Phmmer found no hits! Try HHpred. Exciting...\n")
                return

#     4.3e-36  125.7   6.1    4.6e-36  125.6   6.1    1.0  1  1smm_A    resolution: 1.36 experiment: XRAY release_date: 2004-03-16 [ 348741 : ALL ]
            if CAPTURE:
                hit = line.split()
                if len(hit) >= 9: 
                    if "-----" not in line.split()[0]:
                        if "ECOD" not in DB:
                            if "AFDB" in DB:
                                hitName=hit[8].split(":")[1] + "_PHR"
                            elif "AFCCP4" in DB:
                                hitName=hit[8].split("_")[0].replace("-model","") + "_PHR"
                            else:
                                hitName=hit[8] + "_PHR"
                            TEMPresultsDict[hitName] = PHHit()
                            if "AFDB" in DB:
                                TEMPresultsDict[hitName].afdbName = hit[8].split(":")[1]
                                TEMPresultsDict[hitName].chainID = "A"
                                TEMPresultsDict[hitName].expdta = "AFDB"
                                #TEMPresultsDict[hitName].resolution = float(hit[10])
                                #TEMPresultsDict[hitName].releaseDate = hit[14]
                            elif "AFCCP4" in DB:
                                TEMPresultsDict[hitName].afdbName = hit[8].split("_")[0].replace("-model","")
                                TEMPresultsDict[hitName].chainID = "A"
                                TEMPresultsDict[hitName].expdta = "AFDB"
                            elif "PDBCCP4" in DB:
                                TEMPresultsDict[hitName].afdbName = hit[8][0:4]
                                if len(hit[8]) >= 6:
                                    TEMPresultsDict[hitName].chainID = hit[8][5:]
                                else:
                                    TEMPresultsDict[hitName].chainID = "A"
                                tempRange=hit[20].replace("['", "").replace("']","").split("-")
                                if len(tempRange) == 2:
                                    TEMPresultsDict[hitName].modelResStart= int(tempRange[-2])
                                    TEMPresultsDict[hitName].modelResEnd  = int(tempRange[-1])
                                elif len(tempRange) == 3:
                                    TEMPresultsDict[hitName].modelResStart= -(int(tempRange[-2]))
                                    TEMPresultsDict[hitName].modelResEnd  = int(tempRange[-1])
                                elif len(tempRange) == 4:
                                    TEMPresultsDict[hitName].modelResStart= -(int(tempRange[-3]))
                                    TEMPresultsDict[hitName].modelResEnd  = -(int(tempRange[-1]))
                            else:
                                TEMPresultsDict[hitName].afdbName = hit[8][0:4]
                                if len(hit[8]) >= 6:
                                    TEMPresultsDict[hitName].chainID = hit[8][5:]
                                else:
                                    TEMPresultsDict[hitName].chainID = "A"
                                #tempRange=hit[20].replace("['", "").replace("']","").split("-")
                                #if len(tempRange) == 2:
                                #    TEMPresultsDict[hitName].modelResStart= int(tempRange[-2])
                                #    TEMPresultsDict[hitName].modelResEnd  = int(tempRange[-1])
                                #elif len(tempRange) == 3:
                                #    TEMPresultsDict[hitName].modelResStart= -(int(tempRange[-2]))
                                #    TEMPresultsDict[hitName].modelResEnd  = int(tempRange[-1])
                                #elif len(tempRange) == 4:
                                #    TEMPresultsDict[hitName].modelResStart= -(int(tempRange[-3]))
                                #    TEMPresultsDict[hitName].modelResEnd  = -(int(tempRange[-1]))
                        else:
                            hitName=hit[8].split("|")[1][1:5] + "_" + hit[8].split("|")[1][5:] + "_PHR"
                            TEMPresultsDict[hitName] = PHHit()
                            TEMPresultsDict[hitName].afdbName = hit[8].split("|")[1][1:5]
                            TEMPresultsDict[hitName].chainID = hit[8][5:]
                        TEMPresultsDict[hitName].rank = count
                        TEMPresultsDict[hitName].chainName = hitName 
                        TEMPresultsDict[hitName].score = float(hit[1])
                        TEMPresultsDict[hitName].evalue = float(hit[3])
                        TEMPresultsDict[hitName].ndomains = int(hit[7])
                        #self.resultsDict[hit[1]].prob=float(out[35:40])
                        #self.resultsDict[hit[1]].prob=float(hit[len(hit)-9])
                        #self.resultsDict[hit[1]].score=float(hit[len(hit)-6])
                        #self.resultsDict[hit[1]].evalue=float(hit[len(hit)-8])
                        #self.resultsDict[hit[1]].pvalue=float(hit[len(hit)-7])
                        #self.resultsDict[hit[1]].cols=float(hit[len(hit)-4])
                        count = count + 1
                else:
                    CAPTURE = False

            if "E-value" in line and "score" in line and "bias" in line and "Sequence" in line:
                CAPTURE = True

        # Open the Phmmer log file
        if os.name == "nt":
            plog = open(self.logfile, "r", newline="\r\n")
        else:
            plog = open(self.logfile, "r")
        lines = plog.readlines()
        plog.close()

        import copy

        # Loop over the lines in the log file to grab target-hit alignments
        rawHitList=[]
        count = 0
        for line in lines:
            if ">>" in line[:2]:
                # Check the following line to make sure its not outside the threshold
                if "No individual domains that satisfy reporting thresholds" not in lines[count+1]:
                    if "ECOD" not in DB:
                        if "AFDB" in DB:
                            hit = line.split()[1].split(":")[1] + "_PHR"
                        elif "AFCCP4" in DB:
                            hit = line.split()[1].split("_")[0].replace("-model","") + "_PHR"
                        else:
                            hit = line.split()[1] + "_PHR"
                    else:
                        hit = line.split()[1].split("|")[1][1:5] + "_" + line.split()[1].split("|")[1][5:] + "_PHR"
                    # Count the number of domains
                    domCount = 0
                    domScores=dict([])
                    domline = lines[count + 3]
                    while domline.strip() != "":
                        domCount = domCount + 1
                        if hit in TEMPresultsDict:
                            domainID=int(domline.split()[0])
                            domScores[domainID]=float(domline.split()[2])
                        domline = lines[count + 3 + domCount]
    
                    ecodRange=[]
                    # Grab the alignment and the range it covers in the target sequence
                    for s in range(domCount):
                        if hit in TEMPresultsDict.keys():
                            #hitname = "%s%d" % (hit, s + 1)
                            rawHitList.append(hit)
                            rawCount=rawHitList.count(hit)
                            hitname = "%s_%d_%s" % (hit.split("_")[0], rawCount, hit.split("_")[1])
                            self.resultsList.append(hitname)
    
                            self.resultsDict[hitname] = copy.deepcopy(TEMPresultsDict[hit])
                            self.resultsDict[hitname].domainID = s + 1
                            self.resultsDict[hitname].targetAlignment = (lines[count + 6 + domCount + (6 * s)]).split()[-2].upper()
                            self.resultsDict[hitname].alignment = (lines[count + 8 + domCount + (6 * s)]).split()[-2].upper()
                            self.resultsDict[hitname].score = domScores[self.resultsDict[hitname].domainID]
                            # Get the ranges for the alignment
                            start = (lines[count + 8 + domCount + (6 * s)]).split()[1].upper()
                            end = (lines[count + 8 + domCount + (6 * s)]).split()[-1].upper()
                            if self.resultsDict[hitname].modelResStart is not None:
                                #self.resultsDict[hitname].alnRange = "%d-%d" % (int(start)+self.resultsDict[hitname].modelResStart, int(end)+self.resultsDict[hitname].modelResStart)
                                self.resultsDict[hitname].alnRange = [int(start)+self.resultsDict[hitname].modelResStart, int(end)+self.resultsDict[hitname].modelResStart]
                            else:
                                #self.resultsDict[hitname].alnRange = "%s-%s" % (start, end)
                                self.resultsDict[hitname].alnRange = [int(start.replace("(","")), int(end.replace(")",""))]
                            if "SW" in line.split()[1].split("-")[-1]:
                                self.resultsDict[hitname].modelResRange="[" + line.split("[")[-1]
                            else:
                                self.resultsDict[hitname].modelResRange="['%s-%s']" % (start, end)
                            # If we are using an ECOD database we need to capture all of the ranges presented
                            if "ECOD" in DB:
                                #print(lines[count + 8 + domCount + (6 * s)]).split("|")[-1].split()[0]
                                ecodRange.append((lines[count + 8 + domCount + (6 * s)]).split("|")[-1].split()[0])
                                #ecodRange=(lines[count + 8 + domCount + (6 * s)]).split("|")[-1].split()[0].split(":")[-1]
                                #ecodRange=(lines[count + 8 + domCount + (6 * s)])[0].split("|").split()[-1].split(",")
                                self.resultsDict[hitname].ecodRange=ecodRange
                            startT = (lines[count + 6 + domCount + (6 * s)]).split()[1].upper()
                            endT = (lines[count + 6 + domCount + (6 * s)]).split()[-1].upper()
                            self.resultsDict[hitname].tarRange = [int(startT), int(endT)]
                            self.resultsDict[hitname].tarExtent = (int(endT) - int(startT))
                            self.resultsDict[hitname].tarMidpoint = ((float(endT) - float(startT)) / 2.0) + float(startT)
    
                            simpSID = simpleSeqID.simpleSeqID()
                            local, overall = simpSID.getPercent(self.resultsDict[hitname].alignment,
                                                                self.resultsDict[hitname].targetAlignment, targetSequence)
    
                            self.resultsDict[hitname].localSEQID = local
                            self.resultsDict[hitname].overallSEQID = overall
                            gr = MRBUMP_utils.getPDBres()
                            if self.resultsDict[hitname].expdta != "AFDB":
                                self.resultsDict[hitname].resolution, self.resultsDict[hitname].expdta, self.resultsDict[hitname].releaseDate \
                                    =  gr.getResolution(pdbCODE=self.resultsDict[hitname].afdbName, PDBLOCAL=PDBLOCAL, seqMetaDB=seqMetaDB)
                        
            count = count + 1

        # Figure out the domains for the target that have been matched
        domCount = 1
        self.targetDomainDict[domCount] = Domains()
        self.targetDomainDict[domCount].ID = domCount
        self.targetDomainDict[domCount].midpoint = self.resultsDict[self.resultsList[0]].tarMidpoint
        self.targetDomainDict[domCount].extent = self.resultsDict[self.resultsList[0]].tarExtent
        self.targetDomainDict[domCount].matches.append(self.resultsList[0])
        self.targetDomainDict[domCount].ranges.append(self.resultsDict[self.resultsList[0]].tarRange)

        for hitname in self.resultsList[1:]:
            DOMFOUND = False
            for count in self.targetDomainDict.keys():
                # Has this domain been identified already?
                if self.resultsDict[hitname].tarExtent >= self.targetDomainDict[count].extent - self.extentTolerance and self.resultsDict[hitname].tarExtent <= self.targetDomainDict[count].extent + self.extentTolerance and self.resultsDict[hitname].tarMidpoint >= self.targetDomainDict[count].midpoint - self.midpointTolerance and self.resultsDict[hitname].tarMidpoint <= self.targetDomainDict[count].midpoint + self.midpointTolerance:
                    self.targetDomainDict[count].matches.append(hitname)
                    self.targetDomainDict[count].ranges.append(self.resultsDict[hitname].tarRange)
                    DOMFOUND = True
                    break
            # If we have a new domain set it up
            if DOMFOUND == False:
                domCount = domCount + 1
                self.targetDomainDict[domCount] = Domains()
                self.targetDomainDict[domCount].ID = domCount
                self.targetDomainDict[domCount].midpoint = self.resultsDict[hitname].tarMidpoint
                self.targetDomainDict[domCount].extent = self.resultsDict[hitname].tarExtent
                self.targetDomainDict[domCount].matches.append(hitname)
                self.targetDomainDict[domCount].ranges.append(self.resultsDict[hitname].tarRange)

    def runPhmmer(self,
                  seqin,
                  seqdb,
                  targetSequence,
                  phmmerPickleFile=None,
                  cpu=1,
                  iterations=1,
                  logfile="",
                  alnfile="",
                  debug=False,
                  PDBLOCAL=None,
                  DB=None,
                  seqMetaDB=None,
                  pgap=0.0):
        """ Run Phmmer """

        if logfile == "":
            self.logfile = "phmmer.log"
        else:
            self.logfile = logfile

        if alnfile == "":
            self.alnfile = "phmmerAlignment.log"
        else:
            self.alnfile = alnfile

        self.cpu = cpu
        self.iterations = iterations
        self.seqin = seqin
        self.seqdb = seqdb

        self.termination = False

        # Test that pgap (--pextend) is set correctly
        try:
            pgap = float(pgap)
            if pgap < 0.0 or pgap >= 1.0: 
                sys.stdout.write("pgap (pextend for phmmer) set incorrectly, must be float (0<=x<1.0)")
                sys.exit()
        except ValueError:
            sys.stdout.write("pgap (pextend for phmmer) set incorrectly, must be float (0<=x<1.0)")
            sys.exit()

        # Set the command line
        if pgap != 0.0:
            command_line = self.phmmerEXE + " --pextend %.2lf --cpu %d --notextw --F1 1e-15 --F2 1e-15 -A %s %s %s" % (
                pgap, self.cpu, self.alnfile, self.seqin, self.seqdb)
        else:
            command_line = self.phmmerEXE + " --cpu %d --notextw --F1 1e-15 --F2 1e-15 -A %s %s %s" % (
                self.cpu, self.alnfile, self.seqin, self.seqdb)
        if debug == True:
            sys.stdout.write("Phmmer command line:\n  %s\n" % command_line)
            sys.stdout.write("\n")

        # Launch program
        if os.name == "nt":
            process_args = shlex.split(command_line, posix=False)
            p = subprocess.Popen(
                process_args,
                bufsize=0,
                shell="False",
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)
        else:
            process_args = shlex.split(command_line)
            p = subprocess.Popen(process_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        (child_stdout, child_stdin) = (p.stdout, p.stdin)

        # Write the input
        child_stdin.close()

        # Open the log file for writing
        log = open(self.logfile, "w")
        phmmerALNLog = []

        # Watch the output for successful termination
        out = child_stdout.readline().decode()

        while out:
            phmmerALNLog.append(out)
            if debug == True:
                sys.stdout.write(out)

            if '[ok]' in out:
                self.termination = True

            if "E-value" in out and "score" in out and "bias" in out and "Sequence" in out:
                out = child_stdout.readline().decode()

            log.write(out)
            log.flush()
            out = child_stdout.readline().decode()

        child_stdout.close()

        # Get the alignemnts for each of the hits
        self.getPhmmerAlignments(targetSequence, phmmerALNLog, PDBLOCAL=PDBLOCAL, DB=DB, seqMetaDB=seqMetaDB)

        if phmmerPickleFile is not None:
            pf=open(phmmerPickleFile, "wb")
            pickle.dump(self.resultsDict,pf,pickle.HIGHEST_PROTOCOL)
            pf.close()


if __name__ == "__main__":

    targetSequence = "MGSSHHHHHHSQDPNSSSMAERFDNLVEGLTEERAMAVILADPDSLERPVDKYMAATRLGASNSEESLDVLIQAAELDPEHLFNRITRRKAIDALGRRKSPKALPSLFKALKCSDEAAVINSVEAITKIDAPLTEADHEKLLEALKGEDIQKRAVIQAFCRLGVPGVINSISPLQDDSNPLVAGAARAYMSKVALQPDGLEVLIPQLVDPIAGRRRSAVIDLGDAGDVTRLEALVTAPVSMSLRARSAFQLVDPDKTCQVPEKYAELITQLLQDNPQQLKLRKEWICDIEPTEIENNLQHRDEARQYGGASSLMAMPKAERMILINEIKEKLWSDYVTHYYLTAVVGLQGLEERSDLIRLALAETIPQYTKSRIAAAWGCLRLGLVDQKPLLEELSVSAFWLPLKWTCQRVLKQLS"
    seqin = os.path.join("example.fasta")
    seqdb = os.path.join("mill.fasta")

    # Read in the database
    from mrbump.tools import MRBUMP_utils
    gr = MRBUMP_utils.getPDBres()
    seqMetaDB=gr.readPDBALL()

    ph = phmmer()
    ph.runPhmmer(
        seqin=seqin,
        seqdb=seqdb,
        targetSequence=targetSequence,
        DB="AFDB",
        seqMetaDB=seqMetaDB)

    for x in ph.resultsDict.keys():
        print(x)
        print(ph.resultsDict[x])
