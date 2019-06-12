'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest
import pytest

from mrparse.mr_classify import MrClassifier


@pytest.mark.skip(reason="Need to develop combined classificier tests that run offline")
def test_generate_local():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    topcons_dir = "/opt/MrParse/data/Q13586/topcons"
    classifier = MrClassifier(seqin, topcons_dir=topcons_dir)
    classifier()
    classification = classifier.classification_prediction
    s = '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------MMMMMMMMMMMMMMMMMMMMM-------------CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC---------------C-CC--------CC-CC-----------------------------------------CCCCCCCCCCCCCCC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    assert classification.annotation == s, classification.annotation


if __name__ == '__main__':
    import sys
    pytest.main([__file__] + sys.argv[1:])