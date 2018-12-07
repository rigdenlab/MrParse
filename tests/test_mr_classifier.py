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
    classification = classifier.get_classification(seqin, topcons_dir=topcons_dir)
    s = '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MMMMMMMMMMMMMMMMMMMMM-------------CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC---------------C-CC--------CC-CC-----------------------------------------CCCCCCCCCCCCCCC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    assert classification.annotation == s, classification.annotation
    
def test_generate_remote():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    classification = get_annotation(seqin)
    s = '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MMMMMMMMMMMMMMMMMMMMM-------------CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC---------------C-CC--------CC-CC-----------------------------------------CCCCCCCCCCCCCCC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    assert classification.annotation == s, classification.annotation
