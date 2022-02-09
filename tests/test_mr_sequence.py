#!/usr/bin/env ccp4-python
import set_mrparse_path
import os
from mrparse.mr_sequence import Sequence, merge_multiple_sequences
from Bio import SeqIO
from Bio.Alphabet import generic_protein
import data_constants


def test_2uvo(test_data):
    seq_info = Sequence(seq_file=test_data.x2uvoA_fasta)
    assert seq_info.nresidues == 171
    assert seq_info.sequence == data_constants.TWOUVO_SEQ
    assert abs(seq_info.molecular_weight - 17131) < 0.1
     
 
def test_2uvo_write(test_data):
    seq_info = Sequence(seq_file=test_data.x2uvoA_fasta)
    filename = 'foo.fasta'
    seq_info.write(filename)
    assert os.path.isfile(filename)
    s = SeqIO.read(filename, 'fasta', alphabet=generic_protein)
    assert len(s) == 171
    os.unlink(filename)
    

def test_2uvo_write_description(test_data):
    seq_info = Sequence(seq_file=test_data.x2uvoA_fasta)
    filename = 'foo.fasta'
    description = 'foo'
    seq_info.write(filename, description=description)
    assert os.path.isfile(filename)
    s = SeqIO.read(filename, 'fasta', alphabet=generic_protein)
    assert len(s) == 171
    assert s.description == description
    os.unlink(filename)
    
 
def test_2uvo_seq():
    seq_info = Sequence(sequence=data_constants.TWOUVO_SEQ)
    assert seq_info.nresidues == 171
    assert seq_info.sequence == data_constants.TWOUVO_SEQ
    assert abs(seq_info.molecular_weight - 17131) < 0.1
 
 
def test_2uvo_seq_write():
    seq_info = Sequence(sequence=data_constants.TWOUVO_SEQ)
    filename = 'foo.fasta'
    seq_info.write(filename)
    assert os.path.isfile(filename)
    s = SeqIO.read(filename, 'fasta', alphabet=generic_protein)
    assert len(s) == 171
    os.unlink(filename)


def test_5hxg_seq(test_data):
    seq_info = merge_multiple_sequences(seq_file=test_data.x5hxg_fasta)
    assert seq_info.nresidues == 351
    assert seq_info.sequence == data_constants.FIVEHXG_SEQ


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
