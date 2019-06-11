#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest

from mrparse.mr_sequence import Sequence

TWOUVO_SEQ = 'ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG'
def test_2uvo(test_data):
    seq_info = Sequence(seq_file=test_data.x2uvoA_fasta)
    assert seq_info.nresidues == 171
    assert seq_info.sequence == TWOUVO_SEQ
    assert abs(seq_info.molecular_weight - 17131) < 0.1
    
def test_2uvo_seq(test_data):
    seq_info = Sequence(sequence=TWOUVO_SEQ)
    assert seq_info.nresidues == 171
    assert seq_info.sequence == TWOUVO_SEQ
    assert abs(seq_info.molecular_weight - 17131) < 0.1

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])