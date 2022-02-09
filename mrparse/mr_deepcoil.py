"""
Created on 18 Oct 2018

@author: jmht

"""
import logging
from pathlib import Path
from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation, NULL_ANNOTATION
from mrparse.mr_util import now, is_exe, run_cmd


THRESHOLD_PROBABILITY = 0.6
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
    
    def __init__(self, seq_info, deepcoil_exe):
        self.seq_info = seq_info
        self.deepcoil_exe = deepcoil_exe
        self.prediction = None
        
    def get_prediction(self):
        logger.debug(f"CCPred starting prediction at: {now()}")
        if not is_exe(self.deepcoil_exe):
            raise RuntimeError(f"Cannot find or execute required Deepcoil script: {self.deepcoil_exe}")
        if self.seq_info.nresidues < DEEPCOIL_MIN_RESIDUES or self.seq_info.nresidues > DEEPCOIL_MAX_RESIDUES:
            raise RuntimeError(
                f"Cannot run Deepcoils as sequence length of {self.seq_info.nresidues} "
                f"is outside Deepoil limits of {DEEPCOIL_MIN_RESIDUES} < {DEEPCOIL_MAX_RESIDUES}")
        scores = probabilites_from_sequence(self.seq_info, self.deepcoil_exe)
        ann = SequenceAnnotation()
        ann.source = 'Deepcoil localhost'
        ann.library_add_annotation(CC)
        ann.scores = scores
        ann.annotation = "".join([CC.symbol if float(p) > THRESHOLD_PROBABILITY else NULL_ANNOTATION.symbol
                                  for p in scores])
        logger.debug(f"CCPred finished prediction at: {now()}")
        self.prediction = ann


def probabilites_from_sequence(seq_info, deepcoil_exe):
    output = run_deepcoil(seq_info, deepcoil_exe)
    aa, probabilities = parse_deepcoil(output)
    return probabilities


def run_deepcoil(seq_info, deepcoil_exe):
    """run deepcoil and return the ouptut file
    
    Currently the deepcoil script has no argument to specify the filename and automatically takes it
    from the header of the fasta. As it will mangle the name to ensure a valid filename is produced,
    it makes it hard to know what the file will be called. We therefore write out our own fasta file
    with a set header so we know what the file will be called.
    """
    name = 'deepcoil_input'
    input_fasta = f'{name}.fasta'
    seq_info.write(input_fasta, 'fasta', description=name)
    cmd = [deepcoil_exe,
           '-i',
           input_fasta]
    run_cmd(cmd)
    out_file = f'{name}.out'
    if not Path(out_file).exists:
        logger.debug(f"Could not find named deepcoil output file: {out_file}")
    Path(input_fasta).unlink()
    return out_file


def parse_deepcoil(outfile):
    with open(outfile) as fh:
        tuples = [line.split() for line in fh.readlines()]
    aa, vals = zip(*tuples)
    return aa, vals
