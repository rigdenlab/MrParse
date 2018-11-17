'''
Created on 17 Nov 2018

@author: jmht
'''

import os

TM_SYMBOL = 'M'
CC_SYMBOL = 'C'
HELIX_SYMBOL = 'H'
BSHEET_SYMBOL = 'E'
UNKNOWN_SYMBOL = 'U'


class SequenceChunk(object):
    def __init__(self, start=None, end=None, stype=None, source=None):
        self.start = start
        self.end = end
        self.stype = stype
        self.source = source
    
    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str
    
    
def get_sequence_chunks(sequence, markers=None, source=None):
    assert sequence and markers
    chunks = []
    chunk = None
    for i, s in enumerate(sequence):
        if s in markers:
            if not chunk:
                chunk = SequenceChunk(start=i, stype=s, source=source)
            elif chunk.stype != s:
                chunk.end = i
                chunks.append(chunk)
                chunk = SequenceChunk(start=i, stype=s, source=source)
        else:
            if chunk:
                chunk.end = i
                chunks.append(chunk)
                chunk = None
    if chunk:
        chunk.end = i
        chunks.append(chunk)
    return chunks


def write_fasta(aa, name, maxwidth=80):
    fname = "%s.fasta" % name
    with open(fname, 'w') as w:
        w.write('>%s%s' % (name, os.linesep))
        for chunk in range(0, len(aa), maxwidth):
            w.write(aa[chunk:chunk + maxwidth] + os.linesep)
        w.write(os.linesep)
    return fname

def read_fasta(fasta_file):
    with open(fasta_file) as fh:
        line = fh.readline()
        assert line[0] == '>'
        seq = ""
        while True:
            line = fh.readline()
            if not line:
                break
            line = line.strip()
            if not line or line[0] == '>':
                break
            seq += line
    return seq
                
            
