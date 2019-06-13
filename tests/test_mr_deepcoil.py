#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest

import logging
from mrparse.mr_sequence import Sequence
from mrparse.mr_deepcoil import CCPred

logging.basicConfig(level=logging.DEBUG)


def test_run(test_data):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    cc = CCPred(seq_info)
    cc.get_prediction()
    annotation = cc.prediction
    assert len(annotation) == 171
    assert annotation.scores[0] > 0.01 # first few are always non-zero
    
if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
