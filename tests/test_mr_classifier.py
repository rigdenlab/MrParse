#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest
import logging
import pytest

from mrparse.mr_sequence import Sequence
from mrparse.mr_classify import MrClassifier


logging.basicConfig(level=logging.DEBUG)

def test_generate_null(test_data):
    """This test does nothing as it turns off all of the predictors. It does just however check that the code
    is consistent and there are no syntax errors"""
    seq_info = Sequence(test_data.x2uvoA_fasta)
    classifier = MrClassifier(seq_info=seq_info, do_ss_predictor=False, do_tm_predictor=False, do_cc_predictor=False)
    classifier()
    assert classifier.classification_prediction is None


if __name__ == '__main__':
    import sys
    pytest.main([__file__] + sys.argv[1:])