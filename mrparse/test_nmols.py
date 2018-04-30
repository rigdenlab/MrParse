#!/usr/bin/env ccp4-python

import sys

sys.path.insert(0, '/opt/mrbump-trunk/include/initialisation')
sys.path.insert(0, '/opt/mrbump-trunk/include/file_info')
sys.path.insert(0, '/opt/mrbump-trunk/include/tools')

import Nmol_calc
import MRBUMP_target_info

mtzin = 'toxd.mtz'
seqin = 'toxd.fasta'
logfile = 'matthews.log'

target_info = MRBUMP_target_info.TargetInfo()
target_info.setTargetSeq_Res(seqFile=seqin)
target_info.setTargetWeight()
molecular_weight = target_info.mol_weight

print("TI MOL WEIGHT ", molecular_weight)

Nmol = Nmol_calc.Nmol_calc()
Nmol.calculate_nmol(molecular_weight, mtzin, seqin, logfile)

print("NMOL ", Nmol.nmol_estimate)
print("VM ", Nmol.Vm_estimate)
print("SOLVENT ", Nmol.solvent_estimate)
print("SCORE TOT ", Nmol.score_tot)
