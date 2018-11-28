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

    assert len(annotation) == 685
    assert annotation.annotations[212] == TM.symbol, annotation.annotations[212] 
    assert annotation.annotations[232] == TM.symbol, annotation.annotations[232] 
    assert annotation[212].score == annotation.scores[212]

    assert len(annotation.scores) == 685
    assert annotation.scores[212] > 0.6, annotation.scores[212] 
    assert annotation.scores[232] > 0.6, annotation.scores[232] 


def test_run():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    tc = Topcons()
    prediction = tc.transmembrane_prediction(seqin)
    assert len(prediction) == 685
    assert prediction[212] == 1.0, prediction[212] 
    assert prediction[232] == 1.0, prediction[232]
