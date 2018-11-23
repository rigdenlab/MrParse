'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest

from mrparse.mr_classify import MrClassifier

def test_generate_local():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    topcons_dir = "/opt/MrParse/data/Q13586/topcons"
    classifier = MrClassifier()
    classification = classifier.get_predictions(seqin, topcons_dir=topcons_dir)
    s = '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MMMMMMMMMMMMMMMMMMMMM-------------CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC---------------C-CC--------CC-CC-----------------------------------------CCCCCCCCCCCCCCC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    assert classification == s, classification
