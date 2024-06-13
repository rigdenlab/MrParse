"""
Created on 14 Nov 2018

@author: jmht
"""
import logging
from pathlib import Path
import re
import shutil
import tarfile

from mrparse.mr_annotation import AnnotationSymbol, SequenceAnnotation
from mrparse.mr_util import now, run_cmd

JPRED_SUBMISSION_EMAIL = 'jens.thomas@liverpool.ac.uk'

logger = logging.getLogger(__name__)

HELIX = AnnotationSymbol()
HELIX.symbol = 'H'
HELIX.stype = 'Alpha Helix'
HELIX.name = 'helix'

SHEET = AnnotationSymbol()
SHEET.symbol = 'E'
SHEET.stype = 'Strand'
SHEET.name = 'strand'


class OutOfTimeException(Exception):
    pass


class JPred(object):
    def __init__(self, seq_info=None):
        self.seq_info = seq_info
        self.jpred_script = str(Path(__file__).parent.resolve().joinpath('scripts', 'jpredapi.pl'))
        self.prediction = None
        self.exception = None

    @staticmethod
    def find_results_file(jpred_rundir):
        if not jpred_rundir.exists():
            raise RuntimeError(f"Cannot find directory: {jpred_rundir}")
        out_concise = [str(f) for f in jpred_rundir.iterdir() if str(f).endswith('.concise')][0]
        return out_concise

    @staticmethod
    def parse_jpred_output(out_concise):
        ss_pred = None
        cc_28 = None
        logger.debug(f'Parsing JPRED concise output: {out_concise}')
        with open(out_concise) as f:
            line = f.readline()
            while line:
                prefix = 'Lupas_28:'
                if line.startswith(prefix):
                    line = line.strip().replace(prefix, '')
                    cc_28 = "".join(line.split(","))
                prefix = 'jnetpred:'
                if line.startswith(prefix):
                    line = line.strip().replace(prefix, '')
                    ss_pred = "".join(line.split(","))
                line = f.readline()
        assert ss_pred and cc_28
        return ss_pred, cc_28

    @staticmethod
    def parse_results_output(output):
        """Parse directory path of JPRED
        
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
        mregx = 'Job results archive is now available at: ?(\S+/\S+\.tar\.gz)'
        match = re.search(mregx, output)
        if not match:
            raise RuntimeError(f"Cannot parse directory path from output: {output}")
        dpath = match.group(1)
        return dpath

    @staticmethod
    def parse_status_url(output):
        """Parse jobid and status url for JPRED
        
        Your job will be submitted with the following parameters:
        file: ../data/Q13586.fasta
        format: seq
        skipPDB: on
        email: jens.thomas@liverpool.ac.uk
        name: jens_test_job


        Created JPred job with jobid: jp_H_5vG49
        You can check the status of the job using the following URL: http://www.compbio.dundee.ac.uk/jpred4/cgi-bin/chklog?jp_H_5vG49
        ...or using 'perl jpredapi.pl status jobid=jp_H_5vG49 getResults=yes checkEvery=60 silent' command
        (Check documentation for more details.)

        """
        mregx = 'Created JPred job with jobid: (\S+)\s+You can check the status of the job using the following URL: ?(http?://\S+)'
        match = re.search(mregx, output)
        if not match:
            raise RuntimeError(f"Cannot parse jobid and status_url from output: {output}")
        jobid = match.group(1)
        status_url = match.group(2)
        return jobid, status_url

    @staticmethod
    def create_annotation(annotation):
        ann = SequenceAnnotation()
        ann.source = 'Jpred server'
        ann.annotation = annotation
        ann.library_add_annotation(HELIX)
        ann.library_add_annotation(SHEET)
        ann.scores = [1.0] * len(annotation)
        return ann

    def get_prediction(self, download_tgz=None, jpred_output=None, cleanup=False):
        """Calculate SS using the online JPRED server
        
        Parameters
        ----------
        download_tgz : str
           A results tar.gz file from the JPRED server - FOR RUNNING  UNIT TESTS
        jpred_output : str
           A results file from the JPRED server - FOR RUNNING  UNIT TESTS
        cleanup : bool
           Delete the downloaded/unpacked results
        """
        results_directory = None
        if not (download_tgz or jpred_output):  # for testing
            if not Path(self.seq_info.sequence_file).exists():
                msg = f"Cannot find JPRED sequence file: {self.seq_info.sequence_file}"
                self.exception = msg
                logger.critical(msg)
                raise RuntimeError(msg)
            try:
                download_tgz = self.run_jpred(self.seq_info.sequence_file)
            except Exception as e:
                logger.critical(e)
                self.exception = e
                raise e
        if not jpred_output and not results_directory:
            results_directory = self.unpack_results(download_tgz)
            jpred_output = self.find_results_file(results_directory)
        ss_pred, _ = self.parse_jpred_output(str(jpred_output))
        if cleanup and results_directory.exists():
            self.cleanup(results_directory)
        self.prediction = self.create_annotation(ss_pred)
        logger.debug(f"JPred finished prediction at: {now()}")
        return self.prediction

    def run_jpred(self, seqin):
        logger.debug(f"JPred starting prediction at: {now()}")
        jobid = self.submit_job(seqin)
        download_tgz = self.get_results(jobid)
        return download_tgz

    def submit_job(self, seqin):
        cmd = [self.jpred_script,
               'submit',
               f'file={seqin}',
               'mode=single',
               'format=fasta',
               'name=ccp4_mrparse_submission',
               'skipPDB=on']
        out = run_cmd(cmd)
        jobid, status_url = self.parse_status_url(out)
        logger.info(f"*** Submitted JPRED job with id {jobid} - check its progress here: {status_url}")
        return jobid

    def get_results(self, jobid):
        """Check results and download from the server"""
        cmd = [self.jpred_script,
               'status',
               f'jobid={jobid}',
               'getResults=yes',
               'checkEvery=10']
        out = run_cmd(cmd)
        download_tgz = self.parse_results_output(out)
        download_tgz = Path(download_tgz).resolve()
        logger.debug(f"JPred results downloaded to: {download_tgz}")
        return download_tgz

    def unpack_results(self, download_tgz):
        job_directory = Path(download_tgz).parent
        with tarfile.open(download_tgz, 'r:*') as tf:
            if not tf.getmembers():
                raise RuntimeError(f'Empty archive: {download_tgz}')
            
            import os
            
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tf, path=job_directory)
        logger.debug(f'Extracted jpred files to: {job_directory}')
        return job_directory

    def cleanup(self, results_directory):
        results_directory = Path(results_directory)
        if results_directory.exists():
            logger.debug(f'Removing jpred results directory: {results_directory}')
            shutil.rmtree(str(results_directory))
