#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest

import logging
import pytest
from mrparse.mr_jpred import JPred
from mrparse.mr_sequence import Sequence

logging.basicConfig(level=logging.DEBUG)


def test_parse(test_data):
    p = JPred(jpred_output=test_data.jpred_output).get_prediction()
    assert len(p) == 800
    print(p.annotation)
    assert p.annotation[19:22] == 'EEE', p.annotation[19:22]
    assert p.annotation[612:615]  == 'H-H', p.annotation[612:615]


@pytest.mark.skip(reason="Skipping querying of JPRED server as this is just a unittest.")
def test_submit(test_data):
    seq_info = Sequence(test_data.Q13586_fasta)
    p = JPred(seq_info=seq_info).get_prediction()
    assert len(p) == 685
    assert p.annotation[212] == 'H', p.annotation[212] 
    assert p.annotation[232] == 'H', p.annotation[232] 
    assert p.annotation[390:397] == 'EEEEEEE', p.annotation[390:397]

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
