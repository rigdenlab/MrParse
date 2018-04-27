#!/usr/bin/env ccp4-python

import os
import phaser
import sys
sys.path.insert(0, '/opt/mrbump-trunk/include/file_info')
sys.path.insert(0, '/opt/mrbump-trunk/include/tools')

import MRBUMP_target_info

class Keywords(object):
    """Dummy class fof keywords"""
    def __init__(self):
        self.col_labels = {'F': None, 'SIGF' : None, 'FREER_FLAG' : None}
        self.JOBID = None
        self.AMORE = False
        self.FIXED = False

class Init(object):
    """Dummy init object for testing"""
    def __init__(self, keywords):
        self.hklin = None
        self.init = None
        self.keywords = keywords
        self.ONLYMODELS = False
        self.search_dict = {}
        self.seqin = None

keywords = Keywords()
keywords.JOBID = 'jmht'
keywords.col_labels['F'] =  'FTOXD3'
keywords.col_labels['SIGF'] = 'SIGFTOXD3'
keywords.col_labels['FREER_FLAG'] = 'FreeR_flag'

init = Init(keywords)
init.seqin = 'toxd.fasta'
init.hklin = 'toxd.mtz'
mr_search_dir = "."

# Need to create logs dir for setTargetMTZinfo
logdir = 'logs'
if not os.path.isdir(logdir):
    os.mkdir(logdir)

target_info = MRBUMP_target_info.TargetInfo()
target_info.setMattCoefLogfile(os.path.join(mr_search_dir, logdir, "matt_coef.log"))
retcode = target_info.setTargetInfo(init, mr_search_dir)

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
