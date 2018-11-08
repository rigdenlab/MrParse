'''
Created on 18 Oct 2018

@author: jmht
'''
import os
import math

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


def split_sequence(sequence, chunk_size=500, overlap=100):
    """Split a sequence into chunks of chunk_size that overlap by overlap"""
    
    def _split_sequence(sequence, chunk_size, overlap):
        return [sequence[i:i+chunk_size] for i in range(0, len(sequence)-overlap, chunk_size-overlap)]
    
    if chunk_size < 1:
        raise Exception("chunk size too small")
    if overlap >= chunk_size:
        raise Exception("overlap too large")
    len_seq = len(sequence)
    if len_seq <= chunk_size:
        return [sequence]
    remaining = len_seq % chunk_size
    if remaining == 0:
        chunked = _split_sequence(sequence, chunk_size, overlap)
    else:
        # As the overlap has to be constant we split the list into a set
        # of chunk_size sized chunks, and 2 at the end that will be split
        # into smaller, but equally sized chunks
        nchunks = math.floor(len_seq / chunk_size)
        cut = int(chunk_size * (nchunks - 1))
        head = sequence[0:cut]
        tail = sequence[cut - overlap:] # include the overlap from the previous sequence
        head = _split_sequence(head, chunk_size, overlap)
        tail = _split_sequence(tail, chunk_size, overlap)
        chunked = head + tail
    return chunked
    

# aa = list(range(23))
# print(split_sequence1(aa, chunk_size=4, overlap=2))


# seq = Sequence(fasta=seqin)
# aa = seq.sequences[0]
# num_aa = len(aa)
# if num_aa >= 500:
#     split_sequence(aa)
# 
# run_deepcoil(seqin)

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
