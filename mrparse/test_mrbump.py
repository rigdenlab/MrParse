#!/usr/bin/env ccp4-python

import os
import phaser
import sys
sys.path.insert(0, '/opt/mrbump-trunk/include/building') # for phase_improve
sys.path.insert(0, '/opt/mrbump-trunk/include/cluster') # for phaserEXE
sys.path.insert(0, '/opt/mrbump-trunk/include/ccp4') # for MRBUMP_pdbmerge
sys.path.insert(0, '/opt/mrbump-trunk/include/output') # for WriteLogfile
sys.path.insert(0, '/opt/mrbump-trunk/include/initialisation')
sys.path.insert(0, '/opt/mrbump-trunk/include/file_info')
sys.path.insert(0, '/opt/mrbump-trunk/include/modelling') # for Align_score
sys.path.insert(0, '/opt/mrbump-trunk/include/seq_align') # for Align_score
sys.path.insert(0, '/opt/mrbump-trunk/include/structures') # for Matches
sys.path.insert(0, '/opt/mrbump-trunk/include/tools')

import MRBUMP_target_info
# from mrbump-trunk/include/initialisation
import MRBUMP_keywords
import MRBUMP_initialise

import Matches
import MRBUMP_phmmer


FASTA_DB = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", "pdb70.txt")

class Init(object):
    """Dummy init object for testing"""
    def __init__(self, keywords):
        self.hklin = None
        self.init = None
        self.keywords = keywords
        self.ONLYMODELS = False
        self.search_dict = {}
        self.seqin = None

# Need to create logs dir for setTargetMTZinfo
logdir = 'logs'
if not os.path.isdir(logdir):
    os.mkdir(logdir)

keywords = MRBUMP_keywords.Keywords()
keywords.JOBID = 'jmht'
keywords.col_labels['F'] =  'FTOXD3'
keywords.col_labels['SIGF'] = 'SIGFTOXD3'
keywords.col_labels['FREER_FLAG'] = 'FreeR_flag'

if True:
    # Use dummy object
    init = Init(keywords)
else:
    os.environ['CCP4_BIN'] = os.path.join(os.path.expandvars("$CCP4"),'bin')
    init = MRBUMP_initialise.Initialise()
    init.search_dict = {}

# Set up input data
init.seqin = 'toxd.fasta'
init.hklin = 'toxd.mtz'
init.keywords.JOBID = 'jmht'
init.keywords.col_labels['F'] =  'FTOXD3'
init.keywords.col_labels['SIGF'] = 'SIGFTOXD3'
init.keywords.col_labels['FREER_FLAG'] = 'FreeR_flag'

mr_search_dir = "."
target_info = MRBUMP_target_info.TargetInfo()
target_info.setMattCoefLogfile(os.path.join(mr_search_dir, logdir, "matt_coef.log"))

def test_PHMMER():
    """"see structures/Matches.py"""
    mstat = Matches.Match_struct()
    mstat.FastaDB = FASTA_DB

    retcode = target_info.setTargetInfo(init, mr_search_dir)

    mstat.phmmerLogfile = os.path.join(init.search_dir, "logs", "phmmer.log")
    mstat.phmmerALNfile = os.path.join(init.search_dir, "logs", "phmmerALN.log")
    mstat.phmmerTBLOUT = os.path.join(init.search_dir, "logs", "phmmerTblout.log")
    mstat.phmmerDOMTBLOUT = os.path.join(init.search_dir, "logs", "phmmerDomTblout.log")
    ph = MRBUMP_phmmer.Phmmer()
    ph.runPhmmer(seqin=target_info.seqfile,
                 seqdb=mstat.FastaDB,
                 phmmerTblout=mstat.phmmerTBLOUT,
                 phmmerDomTblout=mstat.phmmerDOMTBLOUT,
                 targetSequence=target_info.sequence,
                 logfile=mstat.phmmerLogfile,
                 alnfile=mstat.phmmerALNfile,
                 debug=True)
    no_of_hits = len(ph.resultsList)
    hitList = ph.resultsList
    mstat.PHresultsDict = ph.resultsDict
    for i in mstat.PHresultsDict.keys():
        sys.stdout.write("%s %s %.2lf %.2lf\n" % (i,
                                                  mstat.PHresultsDict[i].alignment,
                                                  mstat.PHresultsDict[i].localSEQID,
                                                  mstat.PHresultsDict[i].overallSEQID))
        mstat.phmmerALNList.append((i, mstat.PHresultsDict[i].alignment))
        #mstat.makeI2report("%s %s" % (i, mstat.PHresultsDict[i].alignment), mstat.ALN_report)
    source = "PHR"
    with open(os.path.join(init.search_dir, "logs", "phmmer.log"), "r") as f:
        mstat.phmmerLogString = f.read()
    init.search_dict["Sequence_align_list"] = mstat.phmmerALNList

    # Assign the details of the domains in the target based on what has been found by Phmmer
    target_info.targetDomainDict = ph.targetDomainDict
    target_info.noDomains = len(ph.targetDomainDict.keys())

def test_TargetInfo():
    retcode = target_info.setTargetInfo(init, mr_search_dir)

def test_PHMMER2():
   targetSequence="QPRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG"
   seqin     = os.path.join(os.environ["CCP4"], "examples", "toxd", "toxd.seq")
   seqdb     = FASTA_DB
   phmmerTblout    = os.path.join(os.environ["CCP4_SCR"], "phmmerTblout.txt")
   phmmerDomTblout = os.path.join(os.environ["CCP4_SCR"], "phmmerDomTblout.txt")
   phmmerTblout    = "phmmerTblout.txt"
   phmmerDomTblout = "phmmerDomTblout.txt"

   ph = MRBUMP_phmmer.Phmmer()
   ph.runPhmmer(seqin=seqin,
                seqdb=seqdb,
                phmmerTblout=phmmerTblout,
                phmmerDomTblout=phmmerDomTblout,
                targetSequence=targetSequence)

   print ph.resultsList
   print ph.resultsDict["4bd9_B1"].ndomains

   print ph.targetDomainDict[1]


#test_TargetInfo()
test_PHMMER2()
