#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest

import logging
import os
import shutil
import pytest
from mrparse.mr_jpred import JPred
from mrparse.mr_sequence import Sequence

logging.basicConfig(level=logging.DEBUG)


@pytest.mark.skip(reason="Skipping querying of JPRED server as this is just a unittest.")
def test_submit(test_data):
    seq_info = Sequence(test_data.Q13586_fasta)
    p = JPred(seq_info=seq_info).get_prediction()
    assert len(p) == 685
    assert p.annotation[212] == 'H', p.annotation[212] 
    assert p.annotation[232] == 'H', p.annotation[232] 
    assert p.annotation[390:397] == 'EEEEEEE', p.annotation[390:397]


def test_get_prediction_local(test_data):
    p = JPred().get_prediction(jpred_output=test_data.jpred_concise, cleanup=False)
    assert len(p) == 800
    assert p.annotation[19:22] == 'EEE', p.annotation[19:22]
    assert p.annotation[612:615]  == 'H-H', p.annotation[612:615]


def test_upack_and_parse(test_data):
    # Need to setup a directory to run in
    tgz = os.path.basename(test_data.jpred_tgz)
    dir_name = tgz.split('.')[0]
    os.mkdir(dir_name)
    shutil.copy(test_data.jpred_tgz, dir_name)
    download_tgz = os.path.join(dir_name ,tgz)
    p = JPred().get_prediction(download_tgz=download_tgz)
    assert len(p) == 171
    assert p.annotation[22:25] == '-E-', p.annotation[19:22]


def test_parse_results_output():
    results_output = """Your job status will be checked with the following parameters:
JobId: jp_H_5vG49
getResults: yes
checkEvery: 10 [sec]
Thu Nov 29 15:02:01 2018    --->    Job jp_H_5vG49 finished. Results available at the following URL:
http://www.compbio.dundee.ac.uk/jpred4/results/jp_H_5vG49/jp_H_5vG49.results.html


Will attempt to download results now (using 'wget') from:
http://www.compbio.dundee.ac.uk/jpred4/results/jp_H_5vG49/jp_H_5vG49.tar.gz

Job results archive is now available at: jp_H_5vG49/jp_H_5vG49.tar.gz
"""

    dpath = JPred.parse_results_output(results_output)
    assert dpath == 'jp_H_5vG49/jp_H_5vG49.tar.gz'


def test_parse_submission_output():
    submission_output = """Your job will be submitted with the following parameters:
file: ../data/Q13586.fasta
format: seq
skipPDB: on
email: jens.thomas@liverpool.ac.uk
name: jens_test_job


Created JPred job with jobid: jp_H_5vG49
You can check the status of the job using the following URL: http://www.compbio.dundee.ac.uk/jpred4/cgi-bin/chklog?jp_H_5vG49
...or using 'perl jpredapi status jobid=jp_H_5vG49 getResults=yes checkEvery=60 silent' command
(Check documentation for more details.)"""

    jobid, status_url = JPred.parse_status_url(submission_output)
    assert jobid == 'jp_H_5vG49'
    assert status_url == 'http://www.compbio.dundee.ac.uk/jpred4/cgi-bin/chklog?jp_H_5vG49'


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
