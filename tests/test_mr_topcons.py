#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest

import pytest
from mrparse.mr_topcons import TMPred, TM
from mrparse.mr_sequence import Sequence


def test_parse(test_data):
    tc = TMPred(None)
    prediction, scores = tc.parse_topcons_output(test_data.topcons_output)
    annotation = tc.create_annotation(prediction, scores)

    assert len(annotation) == 685
    assert annotation.annotation[212] == TM.symbol, annotation.annotation[212] 
    assert annotation.annotation[232] == TM.symbol, annotation.annotation[232] 
    assert annotation[212].score == annotation.scores[212]

    assert len(annotation.scores) == 685
    assert annotation.scores[212] > 0.6, annotation.scores[212] 
    assert annotation.scores[232] > 0.6, annotation.scores[232] 


@pytest.mark.skip(reason="TMPred server is missing the WSDL description file so can't be used.")
def test_run(test_data):
    seq_info = Sequence(test_data.Q13586_fasta)
    tc = TMPred(seq_info)
    annotation = tc.get_prediction()
    assert len(annotation) == 685
    assert annotation.annotation[212] == TM.symbol, annotation.annotation[212] 
    assert annotation.annotation[232] == TM.symbol, annotation.annotation[232] 
    assert annotation[212].score == annotation.scores[212]
 
    assert len(annotation.scores) == 685
    assert annotation.scores[212] > 0.6, annotation.scores[212] 
    assert annotation.scores[232] > 0.6, annotation.scores[232]


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])