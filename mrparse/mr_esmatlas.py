"""
Created on 04 Nov 2022

@author: rmk65 & hlasimpk
"""

# Simple function to download a predicted model from the ESMAtlas

from Bio import pairwise2
from collections import OrderedDict
import configparser as ConfigParser
from datetime import datetime
import gemmi
import json
import logging
import numpy as np
import os
from pathlib import Path
import requests
from requests.packages.urllib3.exceptions import InsecureRequestWarning
import warnings

from mrparse.mr_alphafold import PdbModelException, calculate_avg_plddt, calculate_sum_plddt, calculate_quality_h_score
from mrparse.mr_alphafold import get_plddt_regions, remove_residues_below_plddt_threshold, convert_plddt_to_bfactor
import mrparse.mr_hit
from mrparse.mr_sequence import Sequence
from simbad.util.pdb_util import PdbStructure
from mrparse.searchDB import phmmer

logger = logging.getLogger(__name__)
warnings.simplefilter('ignore', InsecureRequestWarning)
MODELS_DIR = Path('models')


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
    if not MODELS_DIR.exists():
        MODELS_DIR.mkdir()

    models = OrderedDict()
    for hit in hits.values():
        mlog = ModelData()
        mlog.hit = hit
        hit._homolog = mlog
        mlog.model_url = url = f"https://api.esmatlas.com/fetchPredictedStructure/{hit.pdb_id}.pdb"
        try:
            mlog.pdb_file, mlog.molecular_weight, \
            mlog.avg_plddt, mlog.sum_plddt, mlog.h_score, \
            mlog.date_made, mlog.plddt_regions = prepare_pdb(hit, plddt_cutoff)
        except PdbModelException as e:
            logger.critical(f"Error processing pdb: {e}")
            continue
        models[mlog.name] = mlog
    return models


def download_model(mgnify_id):
    url = f"https://api.esmatlas.com/fetchPredictedStructure/{mgnify_id}.pdb"
    query = requests.get(url, verify=False)
    return query.text

    
def prepare_pdb(hit, plddt_cutoff):
    """
    Download pdb or take file from cache
    trucate to required residues
    calculate the MW
    """

    logger.info("Retrieving and preparing model: %s" % hit.name)

    if hit.model_url is not None:
        pdb_name = os.path.basename(hit.model_url)
    else:
        pdb_name = f"{hit.name.split('_')[0]}"
    pdb_struct = PdbStructure()
    try:
        pdb_string = download_model(pdb_name)
        pdb_struct.structure = gemmi.read_pdb_string(pdb_string)
        pdb_struct.assert_structure()
        pdb_struct.structure.setup_entities()
        if hit.data_created is not None:
            date_made = hit.data_created
        else:
            for line in pdb_string.split('\n'):
                if line.startswith('HEADER'):
                    date_made = line.split()[-1]
    except RuntimeError:
        # SIMBAD currently raises an empty RuntimeError for download problems.
        raise PdbModelException(f"Error downloading PDB file for: {hit.pdb_id}")
    except AssertionError:
        raise PdbModelException(f"Error downloading PDB file for: {hit.pdb_id}")

    pdb_file = MODELS_DIR.joinpath(pdb_name)
    pdb_struct.save(str(pdb_file))

    seqid_range = range(hit.hit_start, hit.hit_stop + 1)

    pdb_struct.structure = convert_lddt_to_plddt(pdb_struct.structure)
    avg_plddt = calculate_avg_plddt(pdb_struct.structure)
    sum_plddt = calculate_sum_plddt(pdb_struct.structure)
    h_score = calculate_quality_h_score(pdb_struct.structure)
    plddt_regions = get_plddt_regions(pdb_struct.structure, hit.seq_ali)

    # Remove residues below threshold
    if plddt_cutoff is not None:
        pdb_struct.structure = remove_residues_below_plddt_threshold(pdb_struct.structure, int(plddt_cutoff))

    # Convert plddt to bfactor score
    pdb_struct.structure = convert_plddt_to_bfactor(pdb_struct.structure)

    truncated_pdb_name = f"esmfold_{hit.pdb_id}.pdb"
    truncated_pdb_path = MODELS_DIR.joinpath(truncated_pdb_name)
    pdb_struct.save(str(truncated_pdb_path),
                    remarks=[f"PHASER ENSEMBLE MODEL 1 ID 1"])
    return str(truncated_pdb_path), int(pdb_struct.molecular_weight), avg_plddt, sum_plddt, h_score, date_made, plddt_regions

def make_esm_prediction(sequence=None, fasta=None, plddt_cutoff=None):
    """ Make a prediction from the ESM Atlas """

    logger.info("Creating an ESMFold model prediction using the ESM Atlas API")

    # Read the input sequence. Sequence string takes priority if both are provided
    if sequence is not None:
        esm_sequence = sequence.strip()
    elif fasta is not None:
        sq = open(fasta, "r")
        sqlines = sq.readlines()
        sq.close()

        esm_sequence = ""
        for line in sqlines:
            if ">" not in line:
                esm_sequence += line.strip()
    else:
        logger.critical("A sequence file or a sequence string needs to be provided")
        raise PdbModelException("No sequence provides")

    if len(esm_sequence) > 400:
        logger.critical("Sequence must be less than 400 residues")
        raise PdbModelException("File too big")

    logger.debug("Using database version specified in mrparse.config")
    config_file = Path(os.environ["CCP4"], "share", "mrparse", "data", "mrparse.config")
    config = ConfigParser.SafeConfigParser()
    config.read(str(config_file))
    database_version = dict(config.items("Databases"))['esmfold_version']

    url = f'https://api.esmatlas.com/foldSequence/{database_version}/pdb/'
    esm_pdb = requests.post(url, data=esm_sequence, timeout=100.0)

    pdb_struct = PdbStructure()
    try:
        pdb_struct.structure = gemmi.read_pdb_string(esm_pdb.text)
        date_made = datetime.today().strftime('%d-%b-%y').upper()
    except FileNotFoundError:
        raise PdbModelException("Error downloading PDB file")

    pdb_struct.structure = convert_lddt_to_plddt(pdb_struct.structure)
    avg_plddt = calculate_avg_plddt(pdb_struct.structure)
    sum_plddt = calculate_sum_plddt(pdb_struct.structure)
    h_score = calculate_quality_h_score(pdb_struct.structure)

    seq_ali = range(1, len(esm_sequence) + 1)
    plddt_regions = get_plddt_regions(pdb_struct.structure, seq_ali)

    # Remove residues below threshold
    if plddt_cutoff is not None:
        pdb_struct.structure = remove_residues_below_plddt_threshold(pdb_struct.structure, int(plddt_cutoff))

    # Convert plddt to bfactor score
    pdb_struct.structure = convert_plddt_to_bfactor(pdb_struct.structure)

    truncated_pdb_name = f"esmfold_{database_version}.pdb"
    truncated_pdb_path = MODELS_DIR.joinpath(truncated_pdb_name)
    pdb_struct.save(str(truncated_pdb_path),
                    remarks=[f"PHASER ENSEMBLE MODEL 1 ID 1"])

    return str(truncated_pdb_path), int(pdb_struct.molecular_weight), avg_plddt, sum_plddt, h_score, date_made, plddt_regions


def convert_lddt_to_plddt(struct):
    for chain in struct[0]:
        for residue in chain:
            for atom in residue:
                lddt_value = atom.b_iso
                atom.b_iso = lddt_value * 100
    return struct


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sequence', help="Input sequence")
    parser.add_argument('-f', '--fasta', help='Input fasta format sequence file (single sequence)')
    parser.add_argument('-plddt_cutoff', help="Set plddt threshold", default=70)
    args = parser.parse_args()

    esm_model, _, _, _, _, _, _ = make_esm_prediction(sequence=args.sequence,
                                                      fasta=args.fasta,
                                                      plddt_cutoff=args.plddt_cutoff)

    print(f"ESMFold model output to: {esm_model}")
