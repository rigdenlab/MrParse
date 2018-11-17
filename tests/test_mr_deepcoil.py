'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest

from mrparse.mr_deepcoil import coiled_coil_prediction


seqin = '/opt/MrParse/data/Q13586.fasta'

predicition = coiled_coil_prediction(seqin)

assert predicition[24] == 0.019
assert len(predicition) == 685
