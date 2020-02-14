'''
Created on 18 Oct 2018

@author: jmht
'''
from collections import OrderedDict
import copy
import logging
import os
from simbad.util.pdb_util import PdbStructure


class PdbModelException(Exception):
    pass


PDB_BASE_URL = 'https://www.rcsb.org/structure/'
PDB_DIR = 'pdb_files'
HOMOLOGS_DIR = 'homologs'

logger = logging.getLogger(__name__)


class HomologData(object):
    OBJECT_ATTRIBUTES = ['hit', 'region']
    def __init__(self):
        self.ellg = None
        self.frac_scat = None
        self.molecular_weight = None
        self.ncopies = None
        self.pdb_url = None
        self.pdb_file = None
        self.resolution = None
        self.rmsd = None
        self.total_frac_scat = None
        self.total_frac_scat_known = None
        # OBJECT_ATTRIBUTES - need to remove from static_dict
        self.hit = None # mr_hit.SequenceHit
        self.region = None
 
    @property
    def chain_id(self):
        return self._get_child_attr('hit', 'chain_id')
        
    @property
    def length(self):
        return self._get_child_attr('hit', 'length')
    
    @property
    def name(self):
        return self._get_child_attr('hit', 'name')
    
    @property
    def pdb_id(self):
        return self._get_child_attr('hit', 'pdb_id')
    
    @property
    def range(self):
        return (self.query_start, self.query_stop)

    @property
    def region_id(self):
        return self._get_child_attr('hit', 'region_id')
    
    @property
    def region_index(self):
        return self._get_child_attr('hit', 'region_index')
    
    @property
    def score(self):
        return self._get_child_attr('hit', 'score')
    
    @property
    def seq_ident(self):
        seq_id = self._get_child_attr('hit', 'local_sequence_identity')
        if seq_id:
            return seq_id / 100.0
        return None
    
    @property
    def query_start(self):
        return self._get_child_attr('hit', 'query_start')
    
    @property
    def query_stop(self):
        return self._get_child_attr('hit', 'query_stop')
    
    # miscellaneous properties
    @property
    def static_dict(self):
        """Return a self representation with all properties resolved, suitable for JSON""" 
        d = copy.copy(self.__dict__)
        for k in self.__dict__.keys():
            if k in self.OBJECT_ATTRIBUTES:
                d.pop(k)
        # Get all properties
        for name in dir(self.__class__):
            obj = getattr(self.__class__, name)
            if name == 'static_dict':
                # Skip this property to avoid infinite recursion
                continue
            if isinstance(obj, property):
                val = obj.__get__(self, self.__class__)
                d[name] = val
        # Need to add in properties as these aren't included
        # FIX ONCE UPDATED JS TO HANDLE TWO INTS
        d['range'] = "{}-{}".format(*self.range)
        return d
    
    def _get_child_attr(self, child, attr):
        if hasattr(self, child):
            child = getattr(self, child)
            if hasattr(child, attr):
                return getattr(child, attr)
        return None
    
    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


def homologs_from_hits(hits, pdb_dir=None):
    pdb_dir = pdb_dir or PDB_DIR
    if not os.path.isdir(pdb_dir):
        os.mkdir(pdb_dir)
    if not os.path.isdir(HOMOLOGS_DIR):
        os.mkdir(HOMOLOGS_DIR)
    homologs = OrderedDict()
    for hit in hits.values():
        hlog = HomologData()
        hlog.hit = hit
        hit._homolog = hlog
        hlog.pdb_url = PDB_BASE_URL + hit.pdb_id
        try:
            hlog.pdb_file, hlog.molecular_weight, hlog.resolution = prepare_pdb(hit, pdb_dir)
        except PdbModelException as e:
            logger.critical("Error processing hit pdb %s", e.message)
        homologs[hlog.name] = hlog
    return homologs


def prepare_pdb(hit, pdb_dir):
    """
    Download pdb or take file from cache
    trucate to required residues
    calculate the MW

    """
    from ample.util.pdb_edit import _select_residues # import on demand as import v slow

    pdb_name = "{}.pdb".format(hit.pdb_id)
    pdb_file = os.path.join(pdb_dir, pdb_name)
    pdb_struct = PdbStructure()
    if os.path.isfile(pdb_file):
        pdb_struct = pdb_struct.from_file(pdb_file)
    else:
        try:
            pdb_struct = pdb_struct.from_pdb_code(hit.pdb_id)
        except RuntimeError:
            # SIMBAD currently raises an empty RuntimeError for download problems.
            raise PdbModelException("Error downloading PDB file for: {}".format(hit.pdb_id))
        _write_pdb_file(pdb_struct, pdb_file)
        #pdb_struct.save(pdb_file)
    
    resolution = pdb_struct.pdb_input.resolution()
        
    pdb_struct.standardize()
    pdb_struct.select_chain_by_id(hit.chain_id)
    if len(pdb_struct.hierarchy.models()) == 0:
        raise PdbModelException("Hierarchy has no models for pdb_name %s" % pdb_name)
    
    seqid_range = range(hit.hit_start, hit.hit_stop + 1) 
    _select_residues(pdb_struct.hierarchy, tokeep_idx=seqid_range)    
    truncated_pdb_name = "{}_{}_{}-{}.pdb".format(hit.pdb_id, hit.chain_id, hit.hit_start, hit.hit_stop)
    truncated_pdb_path = os.path.join(HOMOLOGS_DIR, truncated_pdb_name)
    pdb_struct.save(truncated_pdb_path)
    return truncated_pdb_path, float(pdb_struct.molecular_weight), resolution


def _write_pdb_file(pdb_struct, pdb_file):
    """Horrible hack to write out the PDB file in a way that keeps the TITLE and REMARK records.
    We really just want to save the original PDB but there doesn't seem to be an easy way to do that."""
    pdb_str = ""
    for  s in pdb_struct.pdb_input.title_section():
        pdb_str += s + "\n"
    for  s in pdb_struct.pdb_input.remark_section():
        pdb_str += s + "\n"
    pdb_str += pdb_struct.pdb_input.as_pdb_string()
    with open(pdb_file, 'w') as w:
        w.write(pdb_str)


def calculate_ellg(homologs, hkl_info):
    """Run PHASER to calculate the eLLG values and update the homolog data
    
    Sourced from: ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py"""
    
    if not (hkl_info.molecular_weight and hkl_info.predicted_ncopies):
        raise RuntimeError("Cannot calculate eLLGs without molecular_weight and predicted ncopies")
    import phaser
    mrinput = phaser.InputMR_DAT()
    mrinput.setHKLI(hkl_info.hklin)
    mrinput.setMUTE(True)
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
    ellginput.setMUTE(True)
    
    # Should calculate MW without the search model so that the total MW will be correct when we add the search model
    ellginput.addCOMP_PROT_MW_NUM(hkl_info.molecular_weight, hkl_info.predicted_ncopies)
    search_models = []
    for hname, d in homologs.items():
        if d.pdb_file and d.seq_ident:
            ellginput.addENSE_PDB_ID(hname, d.pdb_file, d.seq_ident)
            search_models.append(hname)
        else:
            d.ellg = -1 # Set to -1 so that sorting works properly
            logger.warn("Cannot calculate eLLG for homolog {} due to missing data.".format(hname))
    ellginput.addSEAR_ENSE_OR_ENSE_NUM(search_models, 1)
    runellg = phaser.runMR_ELLG(ellginput)
    
    phaser_log = 'phaser1.log'
    with open(phaser_log, 'w') as w:
        w.write(runellg.summary())
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
                    h.ellg = float(eLLG)
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
