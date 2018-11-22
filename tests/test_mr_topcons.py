'''
Created on 17 Nov 2018

@author: jmht
'''

import set_mrparse_path
import conftest


from mrparse.mr_sequence import TM_SYMBOL
from mrparse.mr_topcons import Topcons


def test_parse():
    tc = Topcons()
    results_dir = '/opt/MrParse/data/Q13586/topcons'
    prediction, _ = tc.parse_topcons_output(results_dir)
    prediction = [1.0 if p == TM_SYMBOL else 0.0 for p in prediction]
    
    assert len(prediction) == 685
    assert prediction[212] == 1.0, prediction[212] 
    assert prediction[232] == 1.0, prediction[232] 


def test_run():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    tc = Topcons()
    prediction = tc.transmembrane_prediction(seqin)
    assert len(prediction) == 685
    assert prediction[212] == 1.0, prediction[212] 
    assert prediction[232] == 1.0, prediction[232]
