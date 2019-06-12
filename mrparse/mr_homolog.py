'''
Created on 18 Oct 2018

@author: jmht
'''
import copy
import logging
import os
from simbad.util.pdb_util import PdbStructure


class ModelDownloadException(Exception):
    pass


PDB_BASE_URL = 'https://www.rcsb.org/structure/'
PDB_DOWNLOAD_DIR = 'pdb_downloads'
HOMOLOGS_DIR = 'homologs'

logger = logging.getLogger(__name__)

class HomologData(object):
    def __init__(self):
        self.eLLG = None
        self.frac_scat = None
        self.length = None
        self.seq_ident = None
        self.molecular_weight = None
        self.name = None
        self.ncopies = None
        self.pdb_url = None
        self.pdb_file = None
        self.region = None
        self.rmsd = None
        self.score = None
        self.seqid_start = None
        self.seqid_stop = None
        self.total_frac_scat = None
        self.total_frac_scat_known = None
        self._sequence_hit = None # mr_hit.SequenceHit
        
    @property
    def range(self):
        return (self.seqid_start, self.seqid_stop)

    @property
    def range_as_str(self):
        return "{}-{}".format(*self.range)

    @property
    def region_id(self):
        return self.region.ID

    def json_dict(self):
        """Return a representation of ourselves in json""" 
        d = copy.copy(self.__dict__)
        to_remove = ['_sequence_hit']
        for k in self.__dict__.keys():
            if k in to_remove:
                d.pop(k)
        # Need to add in properties as these aren't included
        d['range'] = self.range_as_str
        return d
    
    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


def homologs_from_hits(hits):
    if not os.path.isdir(PDB_DOWNLOAD_DIR):
        os.mkdir(PDB_DOWNLOAD_DIR)
    if not os.path.isdir(HOMOLOGS_DIR):
        os.mkdir(HOMOLOGS_DIR)
    homologs = {}
    for hit in hits.values():
        hlog = HomologData()
        hlog._sequence_hit = hit
        hit._homolog = hlog
        hlog.name = hit.name
        hlog.score = hit.score
        hlog.seq_ident = hit.local_sequence_identity / 100.0
        hlog.region = hit.regionId
        hlog.length = hit.length
        hlog.seqid_start = hit.tarStart
        hlog.seqid_stop = hit.tarStop
        
        hlog.pdb_url = PDB_BASE_URL + hit.pdbName
        try:
            hlog.pdb_file, hlog.molecular_weight = prepare_pdb(hit)
        except ModelDownloadException as e:
            logger.critical("Error processing hit pdb %s", e.message)
        homologs[hlog.name] = hlog
    return homologs


def prepare_pdb(hit):
    """
    Download pdb or take file from cache
    trucate to required residues
    calculate the MW

    """
    from ample.util.pdb_edit import _select_residues # import on demand as import v slow

    pdb_name = "{}_{}.pdb".format(hit.pdbName, hit.chainID)
    pdb_file = os.path.join(PDB_DOWNLOAD_DIR, pdb_name)
    pdb_struct = PdbStructure()
    if os.path.isfile(pdb_file):
        pdb_struct.from_file(pdb_file)
    else:
        pdb_struct.from_pdb_code(hit.pdbName)
        pdb_struct.standardize()
        pdb_struct.select_chain_by_id(hit.chainID)
        pdb_struct.save(pdb_file)
    if len(pdb_struct.hierarchy.models()) == 0:
        raise ModelDownloadException("Hierarchy has no models for pdb_name %s" % pdb_name)
    
    seqid_range = range(hit.tarStart, hit.tarStop + 1) 
    _select_residues(pdb_struct.hierarchy, tokeep_idx=seqid_range)    
    truncated_pdb_name = "{}_{}_{}-{}.pdb".format(hit.pdbName, hit.chainID, hit.tarStart, hit.tarStop)
    truncated_pdb_path = os.path.join(HOMOLOGS_DIR, truncated_pdb_name)
    pdb_struct.save(truncated_pdb_path)
    return truncated_pdb_path, float(pdb_struct.molecular_weight)

def calculate_ellg(homologs, hkl_info):
    """Run PHASER to calculate the eLLG values and update the homolog data
    
    Sourced from: ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py"""
    
    if not (hkl_info.molecular_weight and hkl_info.predicted_ncopies):
        raise RuntimeError("Cannot calculate eLLGs without molecular_weight and predicted ncopies")
    import phaser
    mrinput = phaser.InputMR_DAT()
    mrinput.setHKLI(hkl_info.hklin)
    if hkl_info.labels.i and hkl_info.labels.sigi:
        mrinput.setLABI_I_SIGI(hkl_info.labels.i, hkl_info.labels.sigi)
    elif hkl_info.labels.f and hkl_info.labels.sigf:
        mrinput.setLABI_F_SIGF(hkl_info.labels.f, hkl_info.labels.sigf)
    else:
        msg = "No flags for intensities or amplitudes have been provided"
        raise RuntimeError(msg)

    datrun = phaser.runMR_DAT(mrinput)
    if not datrun.Success():
        raise RuntimeError("Failed to initialise PHASER input.")
    ellginput = phaser.InputMR_ELLG()
    ellginput.setSPAC_HALL(datrun.getSpaceGroupHall())
    ellginput.setCELL6(datrun.getUnitCell())
    ellginput.setREFL_DATA(datrun.getDATA())
    # Can't mute or no logfile!
    #ellginput.setMUTE(True)
    
    # Should calculate MW without the search model so that the total MW will be correct when we add the search model
    ellginput.addCOMP_PROT_MW_NUM(hkl_info.molecular_weight, hkl_info.predicted_ncopies)
    search_models = []
    for hname, d in homologs.items():
        if d.pdb_file and d.seq_ident:
            ellginput.addENSE_PDB_ID(hname, d.pdb_file, d.seq_ident)
            search_models.append(hname)
        else:
            d.eLLG = -1 # Set to -1 so that sorting works properly
            logger.warn("Cannot calculate eLLG for homolog {} due to missing data.".format(hname))
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
            # Get ncopies
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
