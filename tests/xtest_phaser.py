#!/usr/bin/env ccp4-python

import os
import phaser
import sys


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

def percent_solvent_from_vm(vm):
    """Calculation taken from PHASER source code
    phaser/run/runCCA.cc"""
    vm = float(vm)
    assert vm > 0.0, 'VM needs to be > 0.0 : {:F}'.format(vm)
    return (1 - 1.23/vm) * 100

def calc_vm(cell, symmetry, resolution, molecular_weight):
    i = phaser.InputCCA()
    i.setCELL6([float(x) for x in cell.split()])
    i.setRESO_HIGH(resolution)  # 'setRESO', 'setRESO_AUTO_HIGH', 'setRESO_AUTO_OFF', 'setRESO_HIGH', 'setRESO_LOW'
    i.setSPAC_NAME(symmetry)  # setSPAC_HALL'
    # Add fixed
    i.addCOMP_PROT_MW_NUM(molecular_weight, 1)
    i.setMUTE(True)
    #
    #i.addCOMP_PROT_MW_NUM(13600, idx + 1)
    r = phaser.runCCA(i)

    vm = None
    prob = None
    z = None
    if not r.Success():
        raise RuntimeError(r.ErrorName() + " ERROR :" + r.ErrorMessage())

    vm = r.getBestVM()
    prob = r.getBestProb()
    z = r.getBestZ()
    print "Cell Content Analysis"
    print "Molecular weight of assembly = ", r.getAssemblyMW()
    print "Best Z value = ", z
    print "Best VM value = ", vm
    print "Probability of Best VM = ", prob
    print "Solvent ", percent_solvent_from_vm(vm)

    return vm, prob, z

def test_calc_vm():
    cell = '73.582 38.733 23.189 90.000 90.000 90.000'
    symmetry = 'P212121'
    resolution = 2.3
    molecular_weight = 7071.090
    vm , prob, z = calc_vm(cell, symmetry, resolution, molecular_weight)
    nplaces = 5
    assert round(vm, nplaces) == round(2.33662160613, nplaces)

#test_calc_vm()


cell = "72.8400   73.3500   74.2600  103.4000  109.2000  107.4000"
symmetry = "P1"
resolution = 1.962
fixed_MW = 6130*4
target_MW = 13600

i = phaser.InputCCA()
i.setCELL6([float(x) for x in cell.split()])
i.setRESO_HIGH(resolution)  # 'setRESO', 'setRESO_AUTO_HIGH', 'setRESO_AUTO_OFF', 'setRESO_HIGH', 'setRESO_LOW'
i.setSPAC_NAME(symmetry)  # setSPAC_HALL'
# Add fixed
i.addCOMP_PROT_MW_NUM(fixed_MW, 1)
i.addCOMP_PROT_MW_NUM(target_MW, 6)
i.setMUTE(False)
#
#i.addCOMP_PROT_MW_NUM(13600, idx + 1)
r = phaser.runCCA(i)

vm = r.getBestVM()
prob = r.getBestProb()
z = r.getBestZ()
print "Cell Content Analysis"
print "Molecular weight of assembly = ", r.getAssemblyMW()
print "Best Z value = ", z
print "Best VM value = ", vm
print "Probability of Best VM = ", prob
print "Solvent ", percent_solvent_from_vm(vm)



sys.exit()
best_prob = 0.0
for i in range(20):
    j = i + 1
    print "I ",j
    molecular_weight = (target_MW * j) + fixed_MW
    vm, prob, z = calc_vm(cell, symmetry, resolution, molecular_weight)
    print "Z ",z

    if prob > best_prob:
        best_prob = prob
        best_nmol = j

print "GOT ",best_prob, best_nmol


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