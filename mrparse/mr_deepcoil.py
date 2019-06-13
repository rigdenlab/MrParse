'''
Created on 18 Oct 2018

@author: jmht

'''
import logging
import os
import numpy as np
from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation, NULL_ANNOTATION
from mrparse.mr_util import now, is_exe, run_cmd


THRESHOLD_PROBABILITY = 0.6
DEEPCOIL_SCRIPT = '/opt/DeepCoil/deepcoil.py'
DEEPCOIL_MIN_RESIDUES = 30
DEEPCOIL_MAX_RESIDUES = 500

CC = AnnotationSymbol()
CC.symbol = 'C'
CC.stype = 'Coiled-coil Helix'
CC.name = 'CC'

logger = logging.getLogger(__name__)


class CCPred(object):
    """Class to run Coiled-coil prediction using Deepcoil: https://github.com/labstructbioinf/DeepCoil
    
    As Deepcoil has a length limit of 30 < 500 residues, anything over 500 we need to split into chunks
    
    Notes:
    The sequence code should be numpified
    """
    
    def __init__(self, seq_info):
        self.seq_info = seq_info
        self.prediction = None
        
    def get_prediction(self):
        logger.debug("CCPred starting prediction at: %s" % now())
        if not is_exe(DEEPCOIL_SCRIPT):
            raise RuntimeError("Cannot find or execute required Deepcoil script: {}".format(DEEPCOIL_SCRIPT))
        if self.seq_info.nresidues < DEEPCOIL_MIN_RESIDUES or self.seq_info.nresidues > DEEPCOIL_MAX_RESIDUES:
            raise RuntimeError(
                "Cannot run Deepcoils as sequence length of {} is outside Deepoil limits of {} < {}".format(
                    self.seq_info.nresidues, DEEPCOIL_MIN_RESIDUES, DEEPCOIL_MAX_RESIDUES))
        scores = probabilites_from_sequence(self.seq_info)
        ann = SequenceAnnotation()
        ann.source = 'Deepcoil localhost'
        ann.library_add_annotation(CC)
        ann.scores = scores
        ann.annotation = "".join([CC.symbol if p > THRESHOLD_PROBABILITY else NULL_ANNOTATION.symbol for p in scores])
        logger.debug("CCPred finished prediction at: %s" % now())
        self.prediction = ann


def probabilites_from_sequence(seq_info):
    output = run_deepcoil(seq_info)
    aa, probabilities = parse_deepcoil(output)
    return probabilities


def run_deepcoil(seq_info):
    """run deepcoil and return the ouptut file
    
    Currently the deepcoil script has no argument to specify the filename and automatically takes it
    from the header of the fasta. As it will mangle the name to ensure a valid filename is produced,
    it makes it hard to know what the file will be called. We therefore write out our own fasta file
    with a set header so we know what the file will be called.
    """
    name = 'deepcoil_input'
    input_fasta = '{}.fasta'.format(name)
    seq_info.write(input_fasta, 'fasta', description=name)
    cmd = [DEEPCOIL_SCRIPT,
           '-i',
           input_fasta]
    run_cmd(cmd)
    out_file = '{}.out'.format(name)
    if not os.path.isfile(out_file):
        logger.debug("Could not find named deepcoil output file: %s", out_file)
    os.unlink(input_fasta)
    return out_file


def parse_deepcoil(outfile):
    with open(outfile) as fh:
        tuples = [line.split() for line in fh.readlines()]
    aa, vals = zip(*tuples)
    return np.array(aa), np.array(vals, dtype=np.float)
