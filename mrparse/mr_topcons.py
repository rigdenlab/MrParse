'''
Created on 16 Nov 2018

@author: jmht
'''
import os
import shutil
import subprocess
import time
import zipfile


from mrparse.mr_sequence import TM_SYMBOL


class OutOfTimeException(Exception):
    pass


class Topcons(object):
    
    def __init__(self):
        self.max_poll_time = 120
        script_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),'../scripts')
        self.topcons_script = os.path.join(script_dir, 'topcons2_wsdl.py')
        pass

    @staticmethod
    def parse_topcons_output(results_dir):
        assert os.path.isdir(results_dir), "Cannot find directory: %s" % results_dir
        results_file = os.path.join(results_dir, 'query.result.txt')
        with open(results_file) as fh:
            line = fh.readline()
            while line:
                if line.startswith('TOPCONS predicted topology:'):
                    prediction = fh.readline().strip()
                if line.startswith('Predicted TOPCONS reliability'):
                    fh.readline()
                    line = fh.readline().strip()
                    scores = []
                    while line:
                        try:
                            seqid, score = line.split()
                        except ValueError:
                            break
                        scores.append((int(seqid), float(score)))
                        line = fh.readline().strip()
                line = fh.readline()
        return prediction, scores
    
    @staticmethod
    def XXcalculate_probabilities(prediction, scores):
        DEFAULT_PROBABILITY = 50.0
        seqid, score = scores.pop(0)
        probabilities = []
        for i, pred in enumerate(prediction):
            if pred == TM_SYMBOL:
                if seqid == i + 1:  # Assum seqids start from 1
                    prob = score
                else:
                    prob = DEFAULT_PROBABILITY
            else:
                prob = 0.0
            probabilities.append(prob)
            if seqid == i + 1:
                try:
                    seqid, score = scores.pop(0)
                except IndexError:
                    # End of list
                    pass
        return probabilities
    
    def run_topcons(self, seqin):
        jobid = self.submit_job(seqin)
        start = time.time()
        while True:
            elapsed_time = time.time() - start
            if elapsed_time > self.max_poll_time:
                raise OutOfTimeException("Exceed maximum runtime of: %d" % self.max_poll_time)
            if self.job_finished(jobid):
                break
            time.sleep(1)
        results_dir = self.get_results(jobid)
        return results_dir
    
    def submit_job(self, seqin):
        cmd = [self.topcons_script, '-m', 'submit', '-seq', seqin]
        out = subprocess.check_output(cmd)
        print("SUBMIT GOT OUT ",out)
        jobid = None
        for line in out.split(os.linesep):
            if line.startswith("You have successfully submitted your job"):
                jobid = line.split('=')[1].strip()
        if jobid is None:
            raise RuntimeError("Error submitting topcons job: %s" % out)
        return jobid
        
    
    def job_finished(self, jobid):
        """
        Your job with jobid rst_rUF74H is finished!
        The result file ./rst_rUF74H.zip has been retrieved for jobid rst_rUF74H
        """
        cmd = [self.topcons_script, '-m', 'get', '-jobid', jobid]
        out = subprocess.check_output(cmd)
        print("CHECK GOT OUT ",out)
        _jobid = None
        for line in out.split(os.linesep):
            if line.endswith("finished!"):
                _jobid = line.split()[4]
                if _jobid != jobid:
                    raise RuntimeError("Error collecting topcons job: %s" % out)
                else:
                    return True
            elif line.endswith("Please check you typing!"):
                raise RuntimeError("Incorrect jobid: %s" % jobid)
        return False

    
    def get_results(self, jobid):
        ziparchive = jobid + '.zip'
        if not zipfile.is_zipfile(ziparchive):
            raise RuntimeError('File is not a valid zip archive: {0}'.format(ziparchive))
        zipf = zipfile.ZipFile(ziparchive)
        if not zipf.infolist():
            raise RuntimeError('Empty zip file: {0}'.format(ziparchive))
        zipf.extractall()
        assert os.path.isdir(jobid)
        return os.path.abspath(jobid)
    
    def cleanup(self, results_dir):
        if os.path.isdir(results_dir):
            shutil.rmtree(results_dir)
    
    def transmembrane_prediction(self, seqin):
#         try:
        results_dir = self.run_topcons(seqin)
#         except Exception as e:
#             print("Error running topons: %s" % e)
#             return None
        prediction, scores = self.parse_topcons_output(results_dir)
        self.cleanup()
        return [1.0 if p == TM_SYMBOL else 0.0 for p in prediction]
