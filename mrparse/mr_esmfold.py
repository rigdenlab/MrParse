"""
Created on 04 Nov 2022

@author: rmk65 & hlasimpk
"""

# Simple function to download a predicted model from the ESMAtlas

from collections import OrderedDict
import configparser as ConfigParser
from datetime import datetime
import gemmi
import logging
import os
from pathlib import Path
import requests

from mrparse.mr_alphafold import PdbModelException, calculate_avg_plddt, calculate_sum_plddt, calculate_quality_h_score
from mrparse.mr_alphafold import get_plddt_regions, remove_residues_below_plddt_threshold, convert_plddt_to_bfactor
from simbad.util.pdb_util import PdbStructure

logger = logging.getLogger(__name__)
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
        self.range = None

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


def models_from_hits(hit, plddt_cutoff):
    if not MODELS_DIR.exists():
        MODELS_DIR.mkdir()

    models = OrderedDict()
    mlog = ModelData()
    mlog.name = 'ESMFold_model'
    mlog.model_id = 'ESMFold_model'
    mlog.range = (1, len(hit))
    mlog.query_start = 1
    mlog.query_stop = len(hit)
    mlog.region_id = 'A'
    mlog.seq_ident = 1
    mlog.length = len(hit)
    mlog.model_url = 'https://api.esmatlas.com/'
    try:
        mlog.pdb_file, mlog.molecular_weight, \
        mlog.avg_plddt, mlog.sum_plddt, mlog.h_score, \
        mlog.date_made, mlog.plddt_regions = get_esm_prediction(sequence=hit, plddt_cutoff=plddt_cutoff)
    except PdbModelException as e:
        logger.critical(f"Error processing pdb: {e}")
        return {}
    models[mlog.name] = mlog
    return models

def get_esm_prediction(sequence=None, fasta=None, plddt_cutoff=None):
    """ Get a prediction from the ESM Atlas """

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

    if len(esm_sequence) > 350:
        logger.critical("Sequence must be less than 350 residues")
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

    esm_model, _, _, _, _, _, _ = get_esm_prediction(sequence=args.sequence,
                                                     fasta=args.fasta,
                                                     plddt_cutoff=args.plddt_cutoff)

    print(f"ESMFold model output to: {esm_model}")
