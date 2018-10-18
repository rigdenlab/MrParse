#!/usr/bin/env ccp4-python

import os
import phaser
from __builtin__ import False

def anisotropy(HKLIN):
    i = phaser.InputMR_DAT()
    i.setHKLI(HKLIN)
    i.setMUTE(True)
    r = phaser.runMR_DAT(i)
    if not r.Success():
        print "Job exit status FAILURE"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()
    i = phaser.InputANO()
    i.setSPAC_HALL(r.getSpaceGroupHall())
    i.setREFL_DATA(r.getREFL_DATA())
    i.setROOT("beta_blip_ano")
    i.setMUTE(True)
    del (r)
    r = phaser.runANO(i)
    if not r.Success():
        print "Job exit status FAILURE"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()
        return False
    
    print "Anisotropy Correction"
    print "Data read from: ", HKLIN
    print "Data output to : ", r.getMtzFile()
    print "Spacegroup Name (Hall symbol) = %s (%s)" % \
      (r.getSpaceGroupName(), r.getSpaceGroupHall())
    print "Unitcell = ", r.getUnitCell()
    print "Principal components = ", r.getEigenBs()
    print "Range of principal components = ", r.getAnisoDeltaB()
    print "Wilson Scale = ", r.getWilsonK()
    print "Wilson B-factor = ", r.getWilsonB()
    hkl = r.getMiller()
    f = r.getF()
    sigf = r.getSIGF()
    f_iso = r.getCorrectedF(True)
    sigf_iso = r.getCorrectedSIGF(True)
    corr = r.getCorrection(True)
    nrefl = min(10, hkl.size())
    print "First ", nrefl, " reflections"
    print "%4s %4s %4s %10s %10s %10s %10s %10s" % \
      ("H","K","L","F","SIGF",r.getLaboutF(),r.getLaboutSIGF(),"Corr\'n")
    for i in range(0, nrefl):
        print "%4d %4d %4d %10.4f %10.4f %10.4f %10.4f %10.4f" % \
          (hkl[i][0],hkl[i][1],hkl[i][2],f[i],sigf[i],f_iso[i],sigf_iso[i],corr[i])


def has_ncs(HKLIN, F=None, SIGF=None, I=None, SIGI=None):
    i = phaser.InputMR_DAT()
    i.setHKLI(HKLIN)
    useI = False
    useF = False
    if F and SIGF:
        i.setLABI_F_SIGF(F, SIGF)
        useF = True
    elif I and SIGI:
        i.setLABI_I_SIGI(I, SIGI)
        useI = True
    else:
        raise RuntimeError("NEED LABELS!")
    i.setMUTE(True)
    r = phaser.runMR_DAT(i)
    if not r.Success():
        print "Job exit status FAILURE 1"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()
        return False
    i = phaser.InputNCS()
    i.setSPAC_HALL(r.getSpaceGroupHall())
    i.setCELL6(r.getUnitCell())
    if useF:
        i.setREFL_F_SIGF(r.getMiller(), r.getF(), r.getSIGF())
    elif useI:
        i.setREFL_I_SIGI(r.getMiller(), r.getIobs(), r.getSigIobs())
    #i.addCOMP_PROT_MW_NUM(28853, 1)
    #i.addCOMP_PROT_MW_NUM(17522, 1)
    i.setMUTE(False)
    del (r)
    r = phaser.runNCS(i)
    if not r.Success():
        print "Job exit status FAILURE 2"
        print r.ErrorName(), "ERROR :", r.ErrorMessage()
        return False
    print "Translational NCS analysis"
    print "Translational NCS present = ", r.hasTNCS()
    if r.hasTNCS():
        print "Translational NCS vecor = ", r.hasTNCS()
    print "Twinned ",r.twinned()
    print "Twinning alpha = ", r.getTwinAlpha()
    print "DIR ",dir(r)
    return True

if True:
    hklin = os.path.expandvars("$CCP4/lib/py2/site-packages/phaser/tutorial/beta_blip_P3221.mtz")
    has_ncs(hklin, F="Fobs", SIGF="Sigma")
else:
    hklin = os.path.join(os.environ["CEXAM"],"data", "1vr7_lr_i.mtz")
    has_ncs(hklin, I="IMEAN", SIGI="SIGIMEAN")

#anisotropy(hklin)