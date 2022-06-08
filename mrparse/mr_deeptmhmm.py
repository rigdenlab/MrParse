"""
Created on 18 Dec 2020

@author: hlasimpk
"""
import logging
from pyjob import cexec
import glob

from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation
from mrparse.mr_util import now

TM_alpha = AnnotationSymbol()
TM_alpha.symbol = 'M'
TM_alpha.name = 'TM'
TM_alpha.stype = 'Transmembrane Alpha Helix'

TM_beta = AnnotationSymbol()
TM_beta.symbol = 'B'
TM_beta.name = 'TM'
TM_beta.stype = 'Transmembrane Beta Sheet'

logger = logging.getLogger(__name__)


class TMPred(object):
    def __init__(self, seq_info, deeptmhmm_exe="biolib"):
        self.seq_info = seq_info
        self.deeptmhmm_exe = deeptmhmm_exe
        self.prediction = None

    @staticmethod
    def prepare_seqin(seqin):
        with open(seqin) as f:
            lines = f.readlines()
        if len(lines[0].split()) < 2:
            lines[0] = lines[0].replace('\n', ' 1\n')
        with open(seqin, "w") as f:
            f.writelines(lines)

    @staticmethod
    def parse_deeptmhmm_output(annotation_file):
        prediction = ""
        with open(annotation_file) as fh:
            line = fh.readline()
            while line:
                if line.startswith('>'):
                    fh.readline()
                    line = fh.readline()
                prediction += line.replace('\n', '')
                line = fh.readline()
        return prediction

    def create_annotation(self, annotation):
        ann = SequenceAnnotation()
        ann.source = 'DeepTMHMM'
        ann.library_add_annotation(TM_alpha)
        ann.library_add_annotation(TM_beta)
        ann.scores = [0 for _ in annotation]
        ann.annotation = annotation
        return ann

    def run_job(self, seqin):
        cmd = [self.deeptmhmm_exe, 'run', 'DTU/DeepTMHMM', '--fasta', seqin]
        cexec(cmd)
        return

    def get_prediction(self):
        logger.debug(f"DeepTMHMM starting prediction at: {now()}")
        self.prepare_seqin(self.seq_info.sequence_file)
        self.run_job(self.seq_info.sequence_file)
        annotation_file = glob.glob("*/*.3line")[0]
        prediction = self.parse_deeptmhmm_output(annotation_file)
        self.prediction = self.create_annotation(prediction)
        logger.debug(f"DeepTMHMM finished prediction at: {now()}")
