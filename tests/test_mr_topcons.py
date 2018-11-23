'''
Created on 17 Nov 2018

@author: jmht
'''

import set_mrparse_path
import conftest


from mrparse.mr_topcons import Topcons, TM


def test_parse():
    tc = Topcons()
    results_dir = '/opt/MrParse/data/Q13586/topcons'
    prediction, scores = tc.parse_topcons_output(results_dir)
    annotation = tc.create_annotation(prediction, scores)

    assert len(annotation.annotation) == 685
    assert annotation.annotation[212] == TM.symbol, annotation.annotation[212] 
    assert annotation.annotation[232] == TM.symbol, annotation.annotation[232] 

    assert len(annotation.probabilties) == 685
    assert annotation.probabilties[212] > 0.6, annotation.probabilties[212] 
    assert annotation.probabilties[232] > 0.6, annotation.probabilties[232] 


def test_run():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    tc = Topcons()
    prediction = tc.transmembrane_prediction(seqin)
    assert len(prediction) == 685
    assert prediction[212] == 1.0, prediction[212] 
    assert prediction[232] == 1.0, prediction[232]
