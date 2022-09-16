"""
Created on 23 Jul 2021

@author: hlasimpk
"""
from collections import OrderedDict
import configparser as ConfigParser
import ftplib
import gemmi
from itertools import groupby
import logging
import numpy as np
from operator import itemgetter
import os
from pathlib import Path
from pkg_resources import parse_version
import requests
from simbad.util.pdb_util import PdbStructure


class PdbModelException(Exception):
    pass


AF_BASE_URL = 'https://alphafold.ebi.ac.uk/entry/'
AF2_DIR = Path('AF2_files')
MODELS_DIR = Path('models')

logger = logging.getLogger(__name__)


class ModelData(object):
    OBJECT_ATTRIBUTES = ['hit', 'region']

    def __init__(self):
        self.avg_plddt = None
        self.sum_plddt = None
        self.date_made = None
        self.molecular_weight = None
        self.model_url = None
        self.pdb_file = None
        self.h_score = None
        self.rmsd = None
        self.hit = None
        self.region = None
        self.plddt_regions = None

    @property
    def length(self):
        return self._get_child_attr('hit', 'length')

    @property
    def name(self):
        return self._get_child_attr('hit', 'name')

    @property
    def model_id(self):
        return self._get_child_attr('hit', 'name')

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
        d = {k: self.__dict__[k] for k in self.__dict__.keys() if k not in self.OBJECT_ATTRIBUTES}
        # Get all properties
        for name in dir(self.__class__):
            obj = getattr(self.__class__, name)
            if name != 'static_dict' and isinstance(obj, property):
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
        out_str = f"Class: {self.__class__}\nData:\n"
        for a in sorted(attrs):
            out_str += f"  {a} : {self.__dict__[a]}\n"
        return out_str


def models_from_hits(hits, plddt_cutoff):
    if not AF2_DIR.exists():
        AF2_DIR.mkdir()
    if not MODELS_DIR.exists():
        MODELS_DIR.mkdir()

    db_ver = get_afdb_version()
    models = OrderedDict()
    for hit in hits.values():
        mlog = ModelData()
        mlog.hit = hit
        hit._homolog = mlog
        mlog.model_url = AF_BASE_URL + (hit.pdb_id.split("-")[0] + "-" + hit.pdb_id.split("-")[1])
        try:
            mlog.pdb_file, mlog.molecular_weight, \
            mlog.avg_plddt, mlog.sum_plddt, mlog.h_score, \
            mlog.date_made, mlog.plddt_regions = prepare_pdb(hit, plddt_cutoff, db_ver)
        except PdbModelException as e:
            logger.critical(f"Error processing pdb: {e}")
            continue
        models[mlog.name] = mlog
    return models


def download_model(pdb_name):
    """Download AlphaFold2 model"""
    url = 'https://alphafold.ebi.ac.uk/files/' + pdb_name
    query = requests.get(url)
    return query.text


def prepare_pdb(hit, plddt_cutoff, database_version):
    """
    Download pdb or take file from cache
    trucate to required residues
    calculate the MW
    """
    #pdb_name = f"AF-{hit.pdb_id}-F1-model_{database_version}.pdb"
    pdb_name = f"{hit.pdb_id}-model_{database_version}.pdb"
    pdb_struct = PdbStructure()
    try:
        pdb_string = download_model(pdb_name)
        pdb_struct.structure = gemmi.read_pdb_string(pdb_string)
        date_made = pdb_string.split('\n')[0].split()[-1]
    except RuntimeError:
        # SIMBAD currently raises an empty RuntimeError for download problems.
        raise PdbModelException(f"Error downloading PDB file for: {hit.pdb_id}")

    pdb_file = AF2_DIR.joinpath(pdb_name)
    pdb_struct.save(str(pdb_file))

    seqid_range = range(hit.hit_start, hit.hit_stop + 1)
    try:
        pdb_struct.select_residues(to_keep_idx=seqid_range)
    except IndexError:
        # SIMBAD occasionally raises an empty IndexError when selecting residues.
        raise PdbModelException(f"Error selecting residues for: {hit.pdb_id}")

    avg_plddt = calculate_avg_plddt(pdb_struct.structure)
    sum_plddt = calculate_sum_plddt(pdb_struct.structure)
    h_score = calculate_quality_h_score(pdb_struct.structure)
    plddt_regions = get_plddt_regions(pdb_struct.structure, hit.seq_ali)

    # Remove residues below threshold
    if plddt_cutoff is not None:
        pdb_struct.structure = remove_residues_below_plddt_threshold(pdb_struct.structure, int(plddt_cutoff))

    # Convert plddt to bfactor score
    pdb_struct.structure = convert_plddt_to_bfactor(pdb_struct.structure)

    truncated_pdb_name = f"{hit.pdb_id}_{database_version}_{hit.hit_start}-{hit.hit_stop}.pdb"
    truncated_pdb_path = MODELS_DIR.joinpath(truncated_pdb_name)
    pdb_struct.save(str(truncated_pdb_path),
                    remarks=[f"PHASER ENSEMBLE MODEL 1 ID {hit.local_sequence_identity}"])
    return str(truncated_pdb_path), int(pdb_struct.molecular_weight), avg_plddt, sum_plddt, h_score, date_made, plddt_regions


def calculate_quality_threshold(struct, plddt_threshold=70):
    res = above_threshold = 0
    for chain in struct[0]:
        for residue in chain:
            res += 1
            if residue[0].b_iso >= plddt_threshold:
                above_threshold += 1
    return (100.0 / res) * above_threshold


def calculate_quality_h_score(struct):
    score = 0
    for i in reversed(range(1, 101)):
        if calculate_quality_threshold(struct, plddt_threshold=i) >= i:
            score = i
            break
    return score


def get_afdb_version():
    """Query the FTP site to find the latest version of the AFDB"""
    try:
        ftp_host = "ftp.ebi.ac.uk"
        ftp_user = "anonymous"
        ftp_pass = ""
        ftp = ftplib.FTP(ftp_host, ftp_user, ftp_pass)
        ftp.cwd('/pub/databases/alphafold/')
        versions = [x for x in ftp.nlst() if x.startswith('v')]
        return max(versions, key=parse_version)
    except ftplib.all_errors as e:
        logger.debug(f"FTP failed with {e}")
        logger.debug("Using database version specified in mrparse.config")
        config_file = Path(os.environ["CCP4"], "share", "mrparse", "data", "mrparse.config")
        config = ConfigParser.SafeConfigParser()
        config.read(str(config_file))
        return dict(config.items("Databases"))['afdb_version']


def get_plddt(struct):
    plddt_values = []
    for chain in struct[0]:
        for residue in chain:
            plddt_values.append(residue[0].b_iso)
    return plddt_values


def get_plddt_regions(struct, seqid_range):
    regions = {}
    plddt_values = get_plddt(struct)
    residues = zip(seqid_range, plddt_values)

    v_low = []
    low = []
    confident = []
    v_high = []

    for i, plddt in residues:
        if plddt < 50:
            v_low.append(i)
        elif 70 > plddt >= 50:
            low.append(i)
        elif 90 > plddt >= 70:
            confident.append(i)
        elif plddt >= 90:
            v_high.append(i)

    regions['v_low'] = _get_regions(v_low)
    regions['low'] = _get_regions(low)
    regions['confident'] = _get_regions(confident)
    regions['v_high'] = _get_regions(v_high)

    return regions


def _get_regions(residues):
    regions = []
    for k, g in groupby(enumerate(residues), lambda x: x[0] - x[1]):
        group = (map(itemgetter(1), g))
        group = list(map(int, group))
        regions.append((group[0], group[-1]))
    return regions


def calculate_avg_plddt(struct):
    plddt_values = get_plddt(struct)
    return sum(plddt_values) / len(plddt_values)


def calculate_sum_plddt(struct):
    plddt_values = get_plddt(struct)
    return sum(plddt_values)


def convert_plddt_to_bfactor(struct):
    for chain in struct[0]:
        for residue in chain:
            for atom in residue:
                plddt_value = atom.b_iso
                atom.b_iso = _convert_plddt_to_bfactor(plddt_value)
    return struct


def _convert_plddt_to_bfactor(plddt):
    lddt = plddt / 100
    if lddt <= 0.5:
        return 657.97  # Same as the b-factor value with an rmsd estimate of 5.0
    rmsd_est = (0.6 / (lddt ** 3))
    bfactor = ((8 * (np.pi ** 2)) / 3.0) * (rmsd_est ** 2)
    return bfactor


def remove_residues_below_plddt_threshold(struct, plddt_cutoff):
    to_remove = []
    for chain in struct[0]:
        for i, residue in enumerate(chain):
            plddt_value = residue[0].b_iso
            if plddt_value < plddt_cutoff:
                to_remove.append(i)

        for i in to_remove[::-1]:
            del chain[i]
    return struct
