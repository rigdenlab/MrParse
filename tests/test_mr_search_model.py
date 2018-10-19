import sys
import pickle
import set_mrparse_path

from mrparse.mr_search_model import DomainFinder, find_hits, get_homologs

import phaser
 
if False:
    dfinder = DomainFinder()
    seqin = '../data/2uvoA.fasta'
    hits = find_hits(seqin)
    domains = dfinder.find_domains_from_hits(hits)
    homologs = get_homologs(hits, domains)
    
homologs = {'1iqb_B_1': {'mw': 9288.35400000003, 'seqid': 0.45, 'pdb': 'pdb_downloads/1iqb_B.pdb'},
            '2n1s_A_1': {'mw': 2984.445999999996, 'seqid': 0.63, 'pdb': 'pdb_downloads/2n1s_A.pdb'},
            '5xdi_A_1': {'mw': 4085.6129999999916, 'seqid': 0.66, 'pdb': 'pdb_downloads/5xdi_A.pdb'},
            '5wuz_A_1': {'mw': 4559.040999999991, 'seqid': 0.58, 'pdb': 'pdb_downloads/5wuz_A.pdb'},
            '1ulm_B_1': {'mw': 9096.976000000037, 'seqid': 0.54, 'pdb': 'pdb_downloads/1ulm_B.pdb'},
            '2x3t_C_1': {'mw': 16895.870800000117, 'seqid': 1.0, 'pdb': 'pdb_downloads/2x3t_C.pdb'},
            '2lb7_A_1': {'mw': 4430.097999999992, 'seqid': 0.56, 'pdb': 'pdb_downloads/2lb7_A.pdb'},
            '2kus_A_1': {'mw': 3453.903999999995, 'seqid': 0.58, 'pdb': 'pdb_downloads/2kus_A.pdb'},
            '1mmc_A_1': {'mw': 3191.8119999999954, 'seqid': 0.6, 'pdb': 'pdb_downloads/1mmc_A.pdb'},
            '1ulk_B_1': {'mw': 13739.197000000082, 'seqid': 0.53, 'pdb': 'pdb_downloads/1ulk_B.pdb'},
            '4wp4_A_1': {'mw': 4724.179999999993, 'seqid': 0.66, 'pdb': 'pdb_downloads/4wp4_A.pdb'}}


# ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py

mrinput = phaser.InputMR_DAT()
HKLIN = '../data/2uvo_pdbredo.mtz'
mrinput.setHKLI(HKLIN)
F = 'FP'
SIGF = 'SIGFP'
mrinput.setLABI_F_SIGF(F, SIGF)

datrun = phaser.runMR_DAT(mrinput)

if not datrun.Success():
    raise RuntimeError("NO SUCCESS")

ellginput = phaser.InputMR_ELLG()
ellginput.setSPAC_HALL(datrun.getSpaceGroupHall())
ellginput.setCELL6(datrun.getUnitCell())
ellginput.setREFL_DATA(datrun.getDATA())


# addCOMP_PROT_MW_NUM(28853,1)
ellginput.addCOMP_PROT_SEQ_NUM('../data/2uvoA.fasta', 4.0)

for hname, hdata in homologs.items(): 
    ellginput.addENSE_PDB_ID(hname, hdata['pdb'], hdata['seqid'])
    ncopies = 1 # need to calculate?
    ellginput.addSEAR_ENSE_NUM(hname, ncopies)

runellg = phaser.runMR_ELLG(ellginput)

stroutput = runellg.logfile()
with open('phaser.log', 'w') as w:
    w.write(stroutput)


# for ens in MRpdbens:
#     ellginput.addENSE_PDB_ID(ens[0], ens[1], float(ens[2]))
#     ellginput.addSEAR_ENSE_NUM(ens[0], 1)
# for seqfname in seqfnames:
#     ellginput.addCOMP_PROT_SEQ_NUM(seqfname, 1)
# #ellginput.setMUTE(debug==False)
# runellg = runMR_ELLG(ellginput)
# stroutput += runellg.logfile()
# ellgs = []
# for ens in MRpdbens:
#     # append expected ELLGs of each component to the list
#     #import code; code.interact(local=locals())
#     ellg = runellg.get_map_ellg_full_resolution()[ens[0]]
#     ellgs.append((ens[0], Roundoff(ellg, 3)))
# ellglst.append(ellgs)
