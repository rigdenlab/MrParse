'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest

from mrparse.mr_jpred import secondary_structure_prediction


def test_parse():
    jpred_dir = '/opt/MrParse/data/Q13586/jpred'
    p = secondary_structure_prediction(jpred_dir)
    assert p.length == 685
    assert p.annotation[212] == 'H', p.annotation[212] 
    assert p.annotation[232] == 'H', p.annotation[232] 
    assert p.annotation[390:397] == 'EEEEEEE', p.annotation[390:397]

