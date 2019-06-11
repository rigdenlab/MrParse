'''
Created on 17 Nov 2018

@author: jmht
'''

import os

from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
from Bio.Alphabet import IUPAC

SUFFIX_TO_TYPE = { 'fasta' : 'fasta' }


class Sequence(object):
    """Class for handling sequence data"""
    def __init__(self, seq_file, sequence_type=None):
        self.sequence_file = seq_file
        self._molecular_weight = None
        if sequence_type is None:
            sequence_type = self._sequence_type_from_filename(seq_file)
            if not sequence_type:
                raise RuntimeError("Cannot determine sequence type from file: {}".format(seq_file))
        self._bio_seq_record = SeqIO.read(seq_file, sequence_type, alphabet=IUPAC.protein)
        self.nresidues = len(self._bio_seq_record)
        self.sequence = str(self._bio_seq_record.seq)
        
    @staticmethod
    def _sequence_type_from_filename(seq_file):
        _, suffix = os.path.splitext(seq_file)
        suffix = suffix.lstrip('.').lower()
        try:
            return SUFFIX_TO_TYPE[suffix]
        except KeyError:
            return None
    
    @property
    def molecular_weight(self):
        if self._molecular_weight is None:
            self._calculate_molecular_weight()
        return self._molecular_weight
    
    def _calculate_molecular_weight(self):
        self._molecular_weight = molecular_weight(self._bio_seq_record.seq, 'protein')
    
#     def write(self, seq_file, sequence_type=None):
#         if sequence_type is None:
#             sequence_type = self._sequence_type_from_filename(seq_file)
#             if not sequence_type:
#                 raise RuntimeError("Cannot determine sequence type from file: {}".format(seq_file))

# def write_fasta(aa, name, maxwidth=80):
#     fname = "%s.fasta" % name
#     with open(fname, 'w') as w:
#         w.write('>%s%s' % (name, os.linesep))
#         for chunk in range(0, len(aa), maxwidth):
#             w.write(aa[chunk:chunk + maxwidth] + os.linesep)
#         w.write(os.linesep)
#     return fname
# 
# def read_fasta(fasta_file):
#     _, ext = os.path.splitext(fasta_file)
#     seq = ""
#     with open(fasta_file) as fh:
#         for i, line in enumerate(fh.readlines()):
#             if i == 0:
#                 assert line[0] == '>'
#                 continue
#             if i == 1 and ext.lower() == '.seq':
#                 continue # seq files have a blank line 
#             line = line.strip()
#             if not line or line[0] == '>':
#                 break
#             seq += line
#     return seq
                
            
