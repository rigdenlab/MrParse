'''
Created on 14 Nov 2018

@author: jmht
'''
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import time

from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation
from mrparse.mr_util import now

PYTHONVERSION = sys.version_info[0]
JPRED_SUBMISSION_EMAIL = 'jens.thomas@liverpool.ac.uk'

logger = logging.getLogger(__name__)
#logger.addHandler(logging.NullHandler()

HELIX = AnnotationSymbol()
HELIX.symbol = 'H'
HELIX.stype = 'Alpha Helix'
HELIX.name = 'helix'

SHEET = AnnotationSymbol()
SHEET.symbol = 'E'
SHEET.stype = 'B-strand'
SHEET.name = 'b-strand'


class OutOfTimeException(Exception):
    pass


class JPred(object):
    def __init__(self, seq_info=None, jpred_output=None):
        self.seq_info = seq_info
        self.jpred_output = jpred_output
        self.poll_time = 1
        self.max_poll_time = 120
        script_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),'../scripts')
        self.jpred_script = os.path.join(script_dir, 'jpredapi')
        self.prediction = None
        self.exception = None

    def parse_jpred_directory(self, jpred_rundir):
        if not os.path.isdir(jpred_rundir):
            raise RuntimeError("Cannot find directory:%s" % jpred_rundir)
        out_concise = [f for f in os.listdir(jpred_rundir) if f.endswith('.concise')][0]
        out_concise = os.path.join(jpred_rundir, out_concise)
        self.parse_jpred_output(out_concise)

    @staticmethod
    def parse_jpred_output(out_concise):
        cc_28 = None
        ss_pred = None
        logger.debug('Parsing JPRED concise output: %s' % out_concise)
        with open(out_concise) as f:
            line = f.readline()
            while line:
                prefix = 'Lupas_28:'
                if line.startswith(prefix):
                    line = line.strip().replace(prefix,'')
                    cc_28 = "".join(line.split(","))
                prefix = 'jnetpred:'
                if line.startswith(prefix):
                    line = line.strip().replace(prefix,'')
                    ss_pred = "".join(line.split(","))
                line = f.readline()
        assert cc_28 and ss_pred
        return ss_pred, cc_28

    @staticmethod
    def create_annotation(annotation):
        ann = SequenceAnnotation()
        ann.source = 'Jpred server'
        ann.annotation = annotation
        ann.library_add_annotation(HELIX)
        ann.library_add_annotation(SHEET)
        ann.scores = [1.0] * len(annotation)
        return ann

    def get_prediction(self):
        if not self.jpred_output: # for testing
            if not os.path.isfile(self.seq_info.sequence_file):
                msg = "Cannot find JPRED sequence file: %s" % self.seq_info.sequence_file
                self.exception = msg
                logger.critical(msg)
                raise RuntimeError(msg)
            try:
                jpred_rundir = self.run_jpred(self.seq_info.sequence_file)
            except Exception as e:
                logger.critical(e)
                self.exception = e
                raise e
            ss_pred, _ = self.parse_jpred_directory(self.jpred_rundir)
        else:
            jpred_rundir = None
            ss_pred, _ = self.parse_jpred_output(self.jpred_output)
        self.cleanup(jpred_rundir)
        self.prediction = self.create_annotation(ss_pred)
        logger.debug("JPred finished prediction at: %s" % now())
        return self.prediction

    def run_jpred(self, seqin):
        logger.debug("JPred starting prediction at: %s" % now())
        jobid = self.submit_job(seqin)
        start = time.time()
        while True:
            elapsed_time = time.time() - start
            if elapsed_time > self.max_poll_time:
                raise OutOfTimeException("Exceed maximum runtime of: %d" % self.max_poll_time)
            results_path = self.get_results(jobid)
            if results_path:
                break
            time.sleep(self.poll_time)
        results_dir = self.unpack_results(results_path)
        return results_dir
    
    def submit_job(self, seqin):
        """
Your job will be submitted with the following parameters:
file: ../data/Q13586.fasta
format: seq
skipPDB: on
email: jens.thomas@liverpool.ac.uk
name: jens_test_job


Created JPred job with jobid: jp_H_5vG49
You can check the status of the job using the following URL: http://www.compbio.dundee.ac.uk/jpred4/cgi-bin/chklog?jp_H_5vG49
...or using 'perl jpredapi status jobid=jp_H_5vG49 getResults=yes checkEvery=60 silent' command
(Check documentation for more details.)

        """
        cmd = [self.jpred_script,
               'submit',
               'file=%s' % seqin,
               'mode=single',
               'format=fasta',
#                'email={}'.format(JPRED_SUBMISSION_EMAIL),
               'name=ccp4_mrparse_submission',
               'skipPDB=on']
        logger.debug("Running cmd: %s", " ".join(cmd))
        optd = {}
        if PYTHONVERSION > 2:
            optd['encoding'] = 'utf-8'
        out = subprocess.check_output(cmd, **optd)
        logger.debug("%s got output: %s", self.jpred_script, out)
        jobid = None
        for line in out.split(os.linesep):
            if line.startswith("Created JPred job with jobid:"):
                jobid = line.split(':')[1].strip()
        if jobid is None:
            raise RuntimeError("Error submitting jpred job: %s" % out)
        return jobid
        
    def get_results(self, jobid):
        """
Your job status will be checked with the following parameters:
JobId: jp_H_5vG49
getResults: yes
checkEvery: 10 [sec]
Thu Nov 29 15:02:01 2018    --->    Job jp_H_5vG49 finished. Results available at the following URL:
http://www.compbio.dundee.ac.uk/jpred4/results/jp_H_5vG49/jp_H_5vG49.results.html


Will attempt to download results now (using 'wget') from:
http://www.compbio.dundee.ac.uk/jpred4/results/jp_H_5vG49/jp_H_5vG49.tar.gz

Job results archive is now available at: jp_H_5vG49/jp_H_5vG49.tar.gz
        """
        cmd = [self.jpred_script,
               'status',
               'jobid=%s' % jobid,
               'getResults=yes',
               'checkEvery=10']
        logger.debug("Running cmd: %s", " ".join(cmd))
        optd = {}
        if PYTHONVERSION > 2:
            optd['encoding'] = 'utf-8'
        out = subprocess.check_output(cmd, **optd)
        logger.debug("Got output: %s", out)
        _jobid = None
        for line in out.split(os.linesep):
            if line.startswith("Job results archive is now available at:"):
                dpath = line.split(':')[1].strip()
                _jobid, _ = dpath.split('/')
                if _jobid != jobid:
                    raise RuntimeError("Error collecting jpred job: %s" % out)
                else:
                    logger.debug("JPred results downloaded to: %s" % dpath)
                    return dpath
        return False

    def unpack_results(self, results_path):
        jobdir, _ = results_path.split('/')
        with tarfile.open(results_path, 'r:*') as tf:
            if not tf.getmembers():
                raise RuntimeError('Empty archive: %s' % results_path)
            tf.extractall(path=jobdir)
        logger.debug('Extracted jpred files to: %s' % jobdir)
        return jobdir

    def cleanup(self, results_dir):
        if results_dir and os.path.isdir(results_dir):
            logger.debug('Removing jpred results directory: %s' % results_dir)
            shutil.rmtree(results_dir)
