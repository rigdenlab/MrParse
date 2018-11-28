'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest

from mrparse.mr_deepcoil import coiled_coil_prediction


def test_run():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    annotation = coiled_coil_prediction(seqin)
    assert annotation.scores[24] == 0.019
    assert len(annotation) == 685
    
