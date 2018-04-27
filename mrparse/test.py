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

if False:
    # Use dummy object
    init = Init(keywords)
else:
    os.environ['CCP4_BIN'] = os.path.join(os.path.expandvars("$CCP4"),'bin')
    init = MRBUMP_initialise.Initialise()
    init.search_dict = {}

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

    # init
    mstat = Matches.Match_struct()
    #mstat.FastaDB = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", "pdb70_ATOMseqs.txt")
    mstat.FastaDB = os.path.join(os.environ["CCP4"], "share", "mrbump", "data", "pdb70.txt")

    retcode = target_info.setTargetInfo(init, mr_search_dir)

    mstat.phmmerLogfile   = os.path.join(init.search_dir, "logs", "phmmer.log")
    mstat.phmmerALNfile   = os.path.join(init.search_dir, "logs", "phmmerALN.log")
    mstat.phmmerTBLOUT    = os.path.join(init.search_dir, "logs", "phmmerTblout.log")
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
    target_info.targetDomainDict=ph.targetDomainDict
    target_info.noDomains=len(ph.targetDomainDict.keys())

def test_TargetInfo():
    retcode = target_info.setTargetInfo(init, mr_search_dir)


test_PHMMER()
sys.exit()







#Anisotropy Correction 	i = InputANO() 	r = runANO(i) 	ResultANO()
#Cell Content Analysis 	i = InputCCA() 	r = runCCA(i) 	ResultCCA()
#Normal Mode Analysis 	i = InputNMA() 	r = runNMA(i) 	ResultNMA()

# i = InputMR_DAT()
# HKLIn = "beta_blip.mtz"
# F = "Fobs"
# SIGF = "Sigma"
# i.setHKLI(HKLIn)
# i.setLABI(F,SIGF)
# i.setMUTE(True)
# r = runMR_DAT(i)



#
# # COMPosition PROTein MW mw NUM num
# icca.addCOMP_PROT_MW_NUM(mw, num)
# # COMPosition PROTein SEQ file NUM num
# i.addCOMP_PROT_FASTA_NUM(file, num) # where sequence is in a fasta file
# i.addCOMP_PROT_SEQ_NUM(sequence, num) # where sequence is a string, one letter code
# # COMPosition ATOM atomtype NUM num
# i.addCOMP_ATOM_NUM(atomtype,num)
# # COMPosition ENSEmble ens FRAC frac
# i.addCOMP_ENSE_FRAC(ens,frac)
#

# cell, symmetry, resolution

# MRBUMP CODE
# mc.setCELL("72.8400   73.3500   74.2600  103.4000  109.2000  107.4000")
# mc.setSYMM("P1")
# mc.setRESO(1.962)
##               target fixed
# mc.runMC_fixed(13600, 6130*4, "matt_coef.log")

#for idx in range(20):

idx = 4
i = phaser.InputCCA()
i.setCELL6((72.8400, 73.3500, 74.2600, 103.4000, 109.2000, 107.4000))
i.setRESO_HIGH(1.962)  # 'setRESO', 'setRESO_AUTO_HIGH', 'setRESO_AUTO_OFF', 'setRESO_HIGH', 'setRESO_LOW'
i.setSPAC_NAME('P1')  # setSPAC_HALL'
# Add fixed
i.addCOMP_PROT_MW_NUM(6130, 4)
i.setMUTE(False)
# Add target

i.addCOMP_PROT_MW_NUM(13600, idx + 1)
r = phaser.runCCA(i)

print "Cell Content Analysis"
print "Molecular weight of assembly = ", r.getAssemblyMW()
print "Best Z value = ", r.getBestZ()
print "Best VM value = ", r.getBestVM()
print "Probability of Best VM = ", r.getBestProb()

if False:
    i = phaser.InputMR_DAT()
    HKLIN = "../betablip/beta_blip.mtz"
    F = "Fobs"
    SIGF = "Sigma"
    i.setHKLI(HKLIN)
    i.setLABI_F_SIGF(F, SIGF)
    i.setMUTE(True)
    r = phaser.runMR_DAT(i)
    if r.Success():
        print(dir(r))
        i = phaser.InputCCA()
        i.setSPAC_HALL(r.getSpaceGroupHall())
        i.setCELL6(r.getUnitCell())
        i.addCOMP_PROT_MW_NUM(28853, 1)
        i.addCOMP_PROT_MW_NUM(17522, 1)
        i.setMUTE(True)
        del(r)
        r = phaser.runCCA(i)
        if r.Success():
            print "Cell Content Analysis"
            print "Molecular weight of assembly = ", r.getAssemblyMW()
            print "Best Z value = ", r.getBestZ()
            print "Best VM value = ", r.getBestVM()
            print "Probability of Best VM = ", r.getBestProb()
        else:
            print "Job exit status FAILURE"
            print r.ErrorName(), "ERROR :", r.ErrorMessage()
    else:
        print "Job exit status FAILURE"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()

# # Run phaser to position native
# phaser_input = phaser.InputMR_AUTO()
# phaser_input.setSPAC_HALL(hall_symbol)
# phaser_input.setCELL6(unit_cell)
# phaser_input.setREFL_F_SIGF(miller_indices, fp_obs, sigfp_obs)
# # Make sure labels same as native_mtz
# phaser_input.setLABI_F_SIGF(f_label,sigf_label)
# phaser_input.setROOT(root_name)
# phaser_input.addENSE_PDB_ID("native", native_pdb, 0.0)
# phaser_input.addCOMP_PROT_MW_NUM(molecular_weight, 1)
# phaser_input.addSEAR_ENSE_NUM("native", 1)
# phaser_input.setMUTE(True)
# phaser_run = phaser.runMR_AUTO(phaser_input)
# del phaser_input
# if not phaser_run.Success():
#     raise RuntimeError("PHASER failed: {0} : {1}".format(phaser_run.ErrorName(), phaser_run.ErrorMessage()))
#
# if phaser_run.foundSolutions():
#     #print "Phaser has found MR solutions"
#     #print "Top LLG = %f" % phaser_run.getTopLLG()
#     #print "Top PDB file = %s" % phaser_run.getTopPdbFile()
#     #print "Top MTZ file = %s" % phaser_run.getTopMtzFile()
#     native_placed_mtz = phaser_run.getTopMtzFile()
# else:
#     raise RuntimeError("PHASER could not place native_pdb")
#
# if cleanup:
#     os.unlink(phaser_run.getTopPdbFile())
