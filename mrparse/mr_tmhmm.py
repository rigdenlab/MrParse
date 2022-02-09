"""
Created on 18 Dec 2020

@author: hlasimpk
"""
import logging
from pathlib import Path
from pyjob import cexec
import glob

from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation
from mrparse.mr_util import now

TM = AnnotationSymbol()
TM.symbol = 'M'
TM.name = 'TM'
TM.stype = 'Transmembrane Helix'

logger = logging.getLogger(__name__)


class TMPred(object):
    def __init__(self, seq_info, tmhmm_exe="tmhmm"):
        self.seq_info = seq_info
        self.tmhmm_exe = tmhmm_exe
        self.prediction = None
        self.tmhmm_model = str(Path(__file__).parent.resolve().joinpath('scripts', 'TMHMM2.0.model'))

    @staticmethod
    def prepare_seqin(seqin):
        with open(seqin) as f:
            lines = f.readlines()
        if len(lines[0].split()) < 2:
            lines[0] = lines[0].replace('\n', ' 1\n')
        with open(seqin, "w") as f:
            f.writelines(lines)

    @staticmethod
    def parse_tmhmm_output(annotation_file, plot_file):
        prediction = ""
        with open(annotation_file) as fh:
            line = fh.readline()
            while line:
                if line.startswith('>'):
                    line = fh.readline()
                prediction += line.replace('\n', '')
                line = fh.readline()

        probabilties = [[], [], []]
        with open(plot_file) as fh:
            line = fh.readline()
            while line:
                if line.startswith('inside membrane outside'):
                    line = fh.readline()
                probabilties[0].append(line.split()[0])
                probabilties[1].append(line.split()[1])
                probabilties[2].append(line.split()[2])
                line = fh.readline()

        probabilties = TMPred.fix_probabilties(prediction, probabilties)
        return prediction, probabilties

    @staticmethod
    def fix_probabilties(prediction, probabilties):
        _probabilities = []
        for i, pred in enumerate(prediction):
            if pred == "i":
                prob = probabilties[0][i]
            elif pred == "M":
                prob = probabilties[1][i]
            else:
                prob = probabilties[2][i]
            _probabilities.append(prob)
        return _probabilities

    def create_annotation(self, annotation, probabilties):
        annotation = annotation.replace('i', 'M')
        ann = SequenceAnnotation()
        ann.source = 'TMHMM'
        ann.library_add_annotation(TM)
        ann.scores = probabilties
        ann.annotation = annotation
        return ann

    def run_job(self, seqin):
        cmd = [self.tmhmm_exe, '-f', seqin, '-m', self.tmhmm_model]
        cexec(cmd)
        return

    def get_prediction(self):
        logger.debug(f"TMHMM starting prediction at: {now()}")
        # If there are no spaces in the fasta header, tmhmm gives an error
        self.prepare_seqin(self.seq_info.sequence_file)
        self.run_job(self.seq_info.sequence_file)
        annotation_file = glob.glob("*.annotation")[0]
        plot_file = glob.glob("*.plot")[0]
        prediction, scores = self.parse_tmhmm_output(annotation_file, plot_file)
        self.prediction = self.create_annotation(prediction, scores)
        logger.debug(f"TMHMM finished prediction at: {now()}")
