"""
Created on 18 Oct 2018

@author: jmht
"""
from collections import OrderedDict
import copy
import logging
import os, sys
import gzip
import shutil
from pathlib import Path
from simbad.util.pdb_util import PdbStructure


class PdbModelException(Exception):
    pass


PDB_BASE_URL = 'https://www.rcsb.org/structure/'
PDB_DIR = Path('pdb_files')
HOMOLOGS_DIR = Path('homologs')

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
        self.hit = None  # mr_hit.SequenceHit
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
        return self.query_start, self.query_stop

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
        indent = "  "
        out_str = f"Class: {self.__class__}\nData:\n"
        for a in sorted(attrs):
            out_str += indent + f"{a} : {self.__dict__[a]}\n"
        return out_str


def homologs_from_hits(hits, pdb_dir=None, pdb_local=None):
    if not HOMOLOGS_DIR.exists():
        HOMOLOGS_DIR.mkdir()
    homologs = OrderedDict()
    for hit in hits.values():
        hlog = HomologData()
        hlog.hit = hit
        hit._homolog = hlog
        hlog.pdb_url = PDB_BASE_URL + hit.pdb_id
        try:
            hlog.pdb_file, hlog.molecular_weight, hlog.resolution = prepare_pdb(hit, pdb_dir, pdb_local)
        except PdbModelException as e:
            logger.critical(f"Error processing hit pdb {e}")
        homologs[hlog.name] = hlog
    return homologs


def prepare_pdb(hit, pdb_dir, pdb_local):
    """
    Download pdb or take file from cache or local PDB mirror
    truncate to required residues
    calculate the MW

    """

    logger.info("Retrieving and preparing model: %s" % hit.name)
  
    localfile=False
    if pdb_local is not None:
        pdb_local_gzfile=os.path.join(pdb_local, hit.pdb_id[1:3].lower(), "pdb" + hit.pdb_id + ".ent.gz")
        if not os.path.isfile(pdb_local_gzfile):
            logger.info("pdb file not found in local mirror (%s)" % pdb_local_gzfile)
            logger.info("attempting to download or take file from directory instead..")
        else:
            localfile=True
            if not PDB_DIR.exists():
                PDB_DIR.mkdir()
            pdb_file = PDB_DIR.joinpath(f"{hit.pdb_id.lower()}.pdb")
            with gzip.open(pdb_local_gzfile, 'rb') as f_in:
                with open(pdb_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    if not localfile:
        if pdb_dir != "None":
            pdb_dir = Path(pdb_dir)
            entry = list(pdb_dir.glob("*/*"))[0]
            prefix, ext = entry.stem[:-4], entry.suffix
            pdb_file = pdb_dir.joinpath(hit.pdb_id[1:3].lower(), f"{prefix}{hit.pdb_id.lower()}.{ext}")
    
        else:
            if not PDB_DIR.exists():
                PDB_DIR.mkdir()
            pdb_file = PDB_DIR.joinpath(f"{hit.pdb_id.lower()}.pdb")

    pdb_struct = PdbStructure()
    if pdb_file.exists():
        pdb_struct = pdb_struct.from_file(str(pdb_file))
    else:
        try:
            pdb_struct = pdb_struct.from_pdb_code(hit.pdb_id)
            pdb_struct.save(pdb_file)
        except RuntimeError:
            # SIMBAD currently raises an empty RuntimeError for download problems.
            raise PdbModelException(f"Error downloading PDB file for: {hit.pdb_id}")

    resolution = pdb_struct.structure.resolution

    pdb_struct.standardize()
    pdb_struct.select_chain_by_id(hit.chain_id)
    try:
        res_ids = [x.seqid.num for x in pdb_struct.structure[0][0]]
    except IndexError:
        raise PdbModelException(f"Error processing PDB file for: {hit.pdb_id}")
    first_res_id = min(res_ids)
    start = hit.hit_start+int(first_res_id) 
    stop = hit.hit_stop+int(first_res_id) 
    to_keep_seqid_range = range(start, stop + 1)

    to_remove = []
    chain = pdb_struct.structure[0][0]
    for i, residue in enumerate(chain):
        if residue.seqid.num not in to_keep_seqid_range:
            to_remove.append(i)
    for i in to_remove[::-1]:
        del chain[i]

    truncated_pdb_name = f"{hit.pdb_id}_{hit.chain_id}_{start}-{stop}.pdb"
    truncated_pdb_path = HOMOLOGS_DIR.joinpath(truncated_pdb_name)
    pdb_struct.save(str(truncated_pdb_path),
                    remarks=[f"PHASER ENSEMBLE MODEL 1 ID {hit.local_sequence_identity}"])
    return str(truncated_pdb_path), int(round(pdb_struct.molecular_weight)), resolution


def calculate_ellg(homologs, hkl_info):
    """Run PHASER to calculate the eLLG values and update the homolog data
    
    Sourced from: ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py"""

    if not (hkl_info.molecular_weight and hkl_info.predicted_ncopies):
        raise RuntimeError("Cannot calculate eLLGs without molecular_weight and predicted ncopies")
    import phaser
    mrinput = phaser.InputMR_DAT()
    mrinput.setHKLI(hkl_info.hklin)
    mrinput.setMUTE(True)

    hkl_info.input_mtz_obj.read_reflections()
    i_column = hkl_info.input_mtz_obj.i
    sigi_column = hkl_info.input_mtz_obj.sigi
    f_column = hkl_info.input_mtz_obj.f
    sigf_column = hkl_info.input_mtz_obj.sigf

    # Check explicitly for negative intensities
    if i_column and sigi_column and all(hkl_info.input_mtz_obj.reflection_file.column_with_label(column).min_value >= 0 for column in (i_column, sigi_column)):   
        mrinput.setLABI_I_SIGI(i_column, sigi_column)
    elif f_column and sigf_column:
        mrinput.setLABI_F_SIGF(f_column, sigf_column)
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
    struct = None
    for hname, d in homologs.items():
        try:
            struct = PdbStructure().from_file(d.pdb_file)
        except:
            continue
        if d.pdb_file and d.seq_ident and struct:
            ellginput.addENSE_PDB_ID(hname, d.pdb_file, d.seq_ident)
            search_models.append(hname)
        else:
            d.ellg = -1  # Set to -1 so that sorting works properly
            logger.warning(f"Cannot calculate eLLG for homolog {hname} due to missing data.")
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
