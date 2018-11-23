'''
Created on 17 Nov 2018

@author: jmht
'''

import os

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
                
            
