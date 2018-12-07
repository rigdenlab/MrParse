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
    _, ext = os.path.splitext(fasta_file)
    seq = ""
    with open(fasta_file) as fh:
        for i, line in enumerate(fh.readlines()):
            if i == 0:
                assert line[0] == '>'
                continue
            if i == 1 and ext.lower() == '.seq':
                continue # seq files have a blank line 
            line = line.strip()
            if not line or line[0] == '>':
                break
            seq += line
    return seq
                
            
