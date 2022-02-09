"""
Created on 16 Nov 2018

@author: jmht
"""

import logging
import os
from pathlib import Path
import shutil
import time
import zipfile

from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation
from mrparse.mr_util import now, run_cmd


class OutOfTimeException(Exception):
    pass


POLL_TIME = 2
MAX_POLL_TIME = 120

TM = AnnotationSymbol()
TM.symbol = 'M'
TM.name = 'TM'
TM.stype = 'Transmembrane Helix'

logger = logging.getLogger(__name__)


class TMPred(object):
    
    def __init__(self, seq_info):
        self.seq_info = seq_info
        self.prediction = None
        self.topcons_script = str(Path(__file__).parent.resolve().joinpath('scripts', 'topcons2_wsdl.py'))
        self.poll_time = POLL_TIME
        self.max_poll_time = MAX_POLL_TIME

    def parse_topcons_directory(self, results_dir):
        assert Path(results_dir).exists(), f"Cannot find directory: {results_dir}"
        results_file = Path(results_dir, 'query.result.txt')
        return self.parse_topcons_output(str(results_file))

    @staticmethod
    def parse_topcons_output(results_file):
        with open(results_file) as fh:
            line = fh.readline()
            while line:
                if line.startswith('TOPCONS predicted topology:'):
                    prediction = fh.readline().strip()
                if line.startswith('Predicted TOPCONS reliability'):
                    fh.readline()
                    line = fh.readline().strip()
                    probabilties = []
                    while line:
                        try:
                            seqid, prob = line.split()
                        except ValueError:
                            break
                        probabilties.append((int(seqid), float(prob)))
                        line = fh.readline().strip()
                line = fh.readline()
        probabilties = TMPred.fix_probabilties(prediction, probabilties)
        return prediction, probabilties
    
    def create_annotation(self, annotation, probabilties):
        ann = SequenceAnnotation()
        ann.source = 'TopCons server'
        ann.library_add_annotation(TM)
        ann.scores = probabilties
        ann.annotation = annotation
        return ann
    
    @staticmethod
    def fix_probabilties(prediction, probabilties):
        DEFAULT_PROBABILITY = 50.0
        seqid, prob = probabilties.pop(0)
        _probabilities = []
        for i, pred in enumerate(prediction):
            if pred == TM.symbol:
                if seqid == i + 1:  # Assum seqids start from 1
                    new_prob = prob
                else:
                    new_prob = DEFAULT_PROBABILITY
            else:
                new_prob = 0.0
            _probabilities.append(new_prob)
            if seqid == i + 1:
                try:
                    seqid, prob = probabilties.pop(0)
                except IndexError:
                    # End of list
                    pass
        # normalise to between 0.0 and 1.0
        probabilities = [p / 100.0 if p > 0.0 else 0.0 for p in _probabilities]
        return probabilities

    @staticmethod
    def cleanup(results_dir):
        if Path(results_dir).exists():
            shutil.rmtree(results_dir)
    
    def run_topcons(self, seqin):
        jobid = self.submit_job(seqin)
        start = time.time()
        while True:
            elapsed_time = time.time() - start
            if elapsed_time > self.max_poll_time:
                raise OutOfTimeException(f"Exceed maximum runtime of: {self.max_poll_time}")
            if self.job_finished(jobid):
                break
            time.sleep(self.poll_time)
        results_dir = self.get_results(jobid)
        return results_dir
    
    def submit_job(self, seqin):
        cmd = [self.topcons_script, '-m', 'submit', '-seq', seqin]
        out = run_cmd(cmd)
        jobid = None
        for line in out.split(os.linesep):
            if line.startswith("You have successfully submitted your job"):
                jobid = line.split('=')[1].strip()
        if jobid is None:
            raise RuntimeError(f"Error submitting topcons job: {out}")
        return jobid
        
    def job_finished(self, jobid):
        """
        Your job with jobid rst_rUF74H is finished!
        The result file ./rst_rUF74H.zip has been retrieved for jobid rst_rUF74H
        """
        cmd = [self.topcons_script, '-m', 'get', '-jobid', jobid]
        out = run_cmd(cmd)
        _jobid = None
        for line in out.split(os.linesep):
            if line.endswith("finished!"):
                _jobid = line.split()[4]
                if _jobid != jobid:
                    raise RuntimeError(f"Error collecting topcons job: {out}")
                else:
                    return True
            elif line.endswith("Please check you typing!"):
                raise RuntimeError(f"Incorrect jobid: {jobid}")
        return False

    
    def get_results(self, jobid):
        ziparchive = jobid + '.zip'
        if not zipfile.is_zipfile(ziparchive):
            raise RuntimeError(f'File is not a valid zip archive: {ziparchive}')
        zipf = zipfile.ZipFile(ziparchive)
        if not zipf.infolist():
            raise RuntimeError(f'Empty zip file: {ziparchive}')
        zipf.extractall()
        assert Path(jobid).exists()
        return Path(jobid).resolve()
    
    def get_prediction(self):
        logger.debug(f"TMPred starting prediction at: {now()}")
        topcons_dir = self.run_topcons(self.seq_info.sequence_file)
        prediction, scores = self.parse_topcons_directory(topcons_dir)
        self.prediction = self.create_annotation(prediction, scores)
        logger.debug(f"TMPred finished prediction at: {now()}")
        self.cleanup(topcons_dir)
