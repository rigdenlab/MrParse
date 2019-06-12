#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest

import os
from mrparse.mr_sequence import Sequence
from Bio import SeqIO
from Bio.Alphabet import IUPAC
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
    s = SeqIO.read(filename, 'fasta', alphabet=IUPAC.protein)
    assert len(s) == 171
    os.unlink(filename)
    

def test_2uvo_seq():
    seq_info = Sequence(sequence=data_constants.TWOUVO_SEQ)
    assert seq_info.nresidues == 171
    assert seq_info.sequence == data_constants.TWOUVO_SEQ
    assert abs(seq_info.molecular_weight - 17131) < 0.1


def test_2uvo_seq_write(test_data):
    seq_info = Sequence(sequence=data_constants.TWOUVO_SEQ)
    filename = 'foo.fasta'
    seq_info.write(filename)
    assert os.path.isfile(filename)
    s = SeqIO.read(filename, 'fasta', alphabet=IUPAC.protein)
    assert len(s) == 171
    os.unlink(filename)


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])