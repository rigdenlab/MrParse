'''
Created on 18 Oct 2018

@author: jmht
'''
import os
import math

import numpy as np
from ample.util.sequence_util import Sequence
from pyjob import cexec
from pyjob.exception import PyJobExecutionError

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
        nchunks = max(1, math.floor(len_seq / chunk_size))
        cut = max(chunk_size, int(chunk_size * (nchunks - 1)))
        head = sequence[0:cut]
        tail = sequence[cut - overlap:] # include the overlap from the previous sequence
        head = _split_sequence(head, chunk_size, overlap)
        tail = _split_sequence(tail, chunk_size, overlap)
        chunked = head + tail
    return chunked


def run_deepcoil(fasta_in):
    with open(fasta_in) as f:
        name = f.readline().strip()[1:]
    cmd = ['/Users/jmht/miniconda2/envs/py36/bin/python',
            '/opt/DeepCoil/deepcoil.py',
            '-i',
            fasta_in]
    try:
        # Need to sent env so that we don't inherit python2 environment
        stdout = cexec(cmd, env={})
    except (OSError, PyJobExecutionError) as e:
        print("Error running command:{}\n\n{}".format(cmd, e))
        raise(e)
    out_file = '{}.out'.format(name)
    if not os.path.isfile(out_file):
        raise RuntimeError("Deepcoil did not produce output file: %s" % out_file)
    return out_file


def parse_deepcoil(outfile):
    with open(outfile) as fh:
        tuples = [line.split() for line in fh.readlines()]
    aa, vals = zip(*tuples)
    return np.array(aa), np.array(vals, dtype=np.float)


def write_fasta(aa, name, maxwidth=80):
    fname = "%s.fasta" % name
    with open(fname, 'w') as w:
        w.write('>%s%s' % (name, os.linesep))
        for chunk in range(0, len(aa), maxwidth):
            w.write(aa[chunk:chunk + maxwidth] + os.linesep)
        w.write(os.linesep)
    return fname


def join_probability_chunks(probabilities_list, aa_list, overlap):
    cc_prob = []
    aa_check = []
    start = 0
    # Skip last as we always take up to the end
    for p, a in zip(probabilities_list[:-1], aa_list[:-1]):
        p = list(p)
        a = list(a)
        inset = int(overlap / 2)
        end = len(p) - inset
        cc_prob += p[start:end]
        aa_check += a[start:end]
        start = inset
    # add last
    cc_prob += list(probabilities_list[-1])[start:]
    aa_check += list(aa_list[-1])[start:]
    return cc_prob


def run_deepcoil_on_chunks(chunks):
    deepcoil_outputs = []
    for i, chunk in enumerate(chunks):
        name = 'foo%d' % i
        fasta_in = write_fasta(chunk, name)
        dout = run_deepcoil(fasta_in)
        deepcoil_outputs.append(dout)
    return deepcoil_outputs


def parse_chunk_probabilites(deepcoil_outputs):
    aa_list = []
    probabilities_list = []
    for outfile in deepcoil_outputs:
        aa, probs = parse_deepcoil(outfile)
        aa_list.append(aa)
        probabilities_list.append(probs)
    return probabilities_list, aa_list


def coiled_coil_probabilities(seqin, chunk_size = 500, overlap = 100):
#     seq_aa = Sequence(fasta=seqin).sequences[0]
#     chunks = split_sequence(seq_aa, chunk_size=chunk_size, overlap=overlap)
#     deepcoil_outputs = run_deepcoil_on_chunks(chunks)
    deepcoil_outputs = ['foo0.out', 'foo1.out']
    probabilities_list, aa_list = parse_chunk_probabilites(deepcoil_outputs)
    return join_probability_chunks(probabilities_list, aa_list, overlap)


def coiled_coil_predicition(seqin):
    probabilities = coiled_coil_probabilities(seqin)
    prediction = [True if p > 0.6 else False for p in probabilities ]
    chunk_indices = minimal_chunks(prediction, min_chunk=6)
    print("GOT INDICES", chunk_indices)
    # Overwrite probabilties array
    return fill_chunks(probabilities, chunk_indices)
    

def minimal_chunks(prediction, min_chunk=6):
    """Chunk a boolean list, so runs of Trues are at least min_chunk long"""
    p = prediction
    outcount = 0
    inchunk = False
    chunks = []
    start = None
    end = None
    for i in range(len(p)):
        chunk = p[i:i+min_chunk]
        if all(chunk):
            if not inchunk:
                inchunk = True
                start = i
        else:
            if inchunk:
                outcount += 1
                if outcount >= min_chunk:
                    end = i
                    chunks.append((start, end))
                    inchunk = False
                    outcount = 0
    return chunks

def fill_chunks(chunks, chunk_indices):
    chunk_indices = iter(chunk_indices)
    inchunk = False
    start, stop = next(chunk_indices)
    for i in range(len(chunks)):
        if inchunk:
            if i == stop:
                inchunk = False
                chunks[i] = False
                try:
                    start, stop = next(chunk_indices)
                except StopIteration:
                    start = stop = -1
                    inchunk = False
            else:
                chunks[i] = True
        else:
            if i == start:
                inchunk = True
                chunks[i] = True
            else:
                chunks[i] = False
    return chunks


# chunks = [0,0,0,1,1,1,1,0,1,1,0,0,0,1,0,0,1,1,1,1,1,0,0,0]
# min_chunks = minimal_chunks(chunks, min_chunk=3)
# print chunks
# print [1 if x else 0 for x in fill_chunks(chunks, min_chunks)]



seqin = '../data/O75410.fasta'
probs = coiled_coil_predicition(seqin)
print("GOT PROBS ",len(probs), probs)
