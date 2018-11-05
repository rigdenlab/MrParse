'''
Created on 18 Oct 2018

@author: jmht
'''
import os

import numpy as np
from ample.util.sequence_util import Sequence
from pyjob import cexec
from pyjob.script import EXE_EXT


def run_deepcoil(seqin, name='foo'):
    seq = Sequence(fasta=seqin)
    seq.headers = ['>{}'.format(name)]
    fasta_in = name + '.fasta'
    seq.write_fasta(fasta_in)
    cmd = ['/Users/jmht/miniconda2/envs/py36/bin/python',
            '/opt/DeepCoil/deepcoil.py',
            '-i',
            fasta_in]
    try:
        cexec(cmd)
    except OSError as e:
        print("Error running command:{}\n{}".format(cmd, e))
    out_file = '{}.out'.format(name)
    return out_file 

def parse_deepcoil(outfile):
    with open(outfile) as fh:
        tuples = [line.split() for line in fh.readlines()]
    aa, vals = zip(*tuples)
    return np.array(aa), np.array(vals, dtype=np.float)
        
        

seqin = '../data/O75410.fasta'

run_deepcoil(seqin)
# seq = Sequence(fasta=seqin)
# 
# d1f = '/opt/DeepCoil/O75410_1_500.out'
# aa1, vals1 = parse_deepcoil(d1f)
# 
# d2f = '/opt/DeepCoil/O75410_400_.out'
# aa2, vals2 = parse_deepcoil(d2f)
# 
# print(aa1[-100:])
# print(aa2[0:100])
# 
# print(vals1[400:500])
# print(vals2[0:100])

# First 5 of 2nd sequence are different (non-zero vs zero) - probably edge effect
# Values differ in last ~25 places
# Overlap by 100 and take 0:-50 of first and 50:- of second sequnce
