'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import os

from simbad.util.pdb_util import PdbStructure

PDB_BASE_URL = 'https://www.rcsb.org/structure/'

logger = logging.getLogger(__name__)

class HomologData(object):
    def __init__(self):
        self.name = None
        self.eLLG = None
        self.frac_scat = None
        self.total_frac_scat = None
        self.total_frac_scat_known = None
        self.rmsd = None
        self.ncopies = None
        self.molecular_weight = None
        self.score = None
        self.seqid = None
        self.pdb_url = None
        self.pdb_file = None
        self.domain = None
        self.range = None

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


def get_homologs(hits, domains):
    pdb_dir = 'pdb_downloads'
    if not os.path.isdir(pdb_dir):
        os.mkdir(pdb_dir)
    homologs = {}
    for domain in domains:
        for match in domain.matches:
            hit = hits[match]
            pdb_name = hit.pdbName + '_' + hit.chainID + '.pdb'
            pdb_file = os.path.join(pdb_dir, pdb_name)
            pdb_struct = PdbStructure()
            if os.path.isfile(pdb_file):
                pdb_struct.from_file(pdb_file)
            else:
                pdb_struct.from_pdb_code(hit.pdbName)
                pdb_struct.standardize()
                pdb_struct.select_chain_by_id(hit.chainID)
                pdb_struct.save(pdb_file)
            hlog = HomologData()
            hlog.name = hit.name
            hlog.score = hits[hit.name].score
            hlog.seqid = hits[hit.name].localSEQID / 100.0
            hlog.domain = domain.ID
            hlog.range = hits[hit.name].tarRange
            homologs[hit.name] = hlog
            hlog.pdb_url = PDB_BASE_URL + hit.pdbName
            if len(pdb_struct.hierarchy.models()) == 0:
                logger.critical("Hierarchy has no models for pdb_name %s" % pdb_name)
            else:
                hlog.pdb_file = pdb_file
                hlog.molecular_weight = float(pdb_struct.molecular_weight)
    return homologs


def ellg_data_from_phaser_log(fpath, homologs):
    with open(fpath) as fh:
        line = fh.readline()
        while line:
            # Get base homolog data
            if line.strip().startswith('eLLG: eLLG of ensemble alone'):
                fh.readline()
                while True:
                    line = fh.readline()
                    if not line.strip():
                        break
                    eLLG, rmsd, frac_scat, name = line.strip().split()
                    h = homologs[name]
                    h.eLLG = float(eLLG)
                    h.rmsd = float(rmsd)
                    h.frac_scat = float(frac_scat)
            # Get ncoopies
            if line.strip().startswith('Number of copies for eLLG target'):
                fh.readline()
                fh.readline()
                while True:
                    line = fh.readline()
                    if not line.strip():
                        break
                    _, _, total_frac_scat_known, total_frac_scat, ncopies, homolog = line.strip().split()
                    h = homologs[homolog]
                    h.total_frac_scat_known = float(total_frac_scat_known)
                    h.total_frac_scat = float(total_frac_scat)
                    h.ncopies = int(ncopies)
            line = fh.readline()
    return homologs


def calculate_ellg(homologs, mtz):
    """Stuff from : ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py"""
    import phaser
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
    # Can't mute or no logfile!
    #ellginput.setMUTE(True)
    
    asu_mw = 72846.44

    # Should calculate MW without the search model so that the total MW will be correct when we add the search model
    ellginput.addCOMP_PROT_MW_NUM(asu_mw, 1)
    search_models = []
    for hname, d in homologs.items():
        ellginput.addENSE_PDB_ID(hname, d.pdb_file, d.seqid)
        search_models.append(hname)
    ellginput.addSEAR_ENSE_OR_ENSE_NUM(search_models, 1)
    
    runellg = phaser.runMR_ELLG(ellginput)
    
    """
    DIR runellg ['ErrorMessage', 'ErrorName', 'Failed', 'Failure', 'Success',  'concise', 'cpu_time', 
    'get_ellg_full_resolution', 'get_ellg_target_nres', 'get_map_chain_ellg', 'get_map_ellg_full_resolution', 
    'get_perfect_data_resolution', 'get_target_resolution', 'get_useful_resolution', 'git_branchname', 'git_commitdatetime', 
    'git_hash', 'git_rev', 'git_shorthash', 'git_tagname', 'git_totalcommits', 'logfile', 'loggraph', 'output_strings', 'process', 
    'run_time', 'setLevel', 'setPackageCCP4', 'setPackagePhenix', 'setPhenixCallback', 'setPhenixPackageCallback', 'set_callback', 
    'set_file_object', 'set_sys_stdout', 'summary', 'svn_revision', 'verbose', 'version_number', 'warnings']
    """
    
    stroutput = runellg.logfile()
    phaser_log = 'phaser1.log'
    with open(phaser_log, 'w') as w:
        w.write(stroutput)
    ellg_data_from_phaser_log(phaser_log, homologs)

