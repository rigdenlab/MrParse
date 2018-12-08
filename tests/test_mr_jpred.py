'''
Created on 17 Nov 2018

@author: jmht
'''
import logging
import set_mrparse_path
import conftest

from mrparse.mr_jpred import JPred

logging.basicConfig(level=logging.DEBUG)


def test_parse():
    jpred_rundir = '/opt/MrParse/data/Q13586/jpred'
    p = JPred(jpred_rundir=jpred_rundir).get_prediction()
    assert len(p) == 685
    assert p.annotation[212] == 'H', p.annotation[212] 
    assert p.annotation[232] == 'H', p.annotation[232] 
    assert p.annotation[390:397] == 'EEEEEEE', p.annotation[390:397]


def test_submit():
    seqin = '/opt/MrParse/data/Q13586.fasta'
    p = JPred(seqin=seqin).get_prediction()
    assert len(p) == 685
    assert p.annotation[212] == 'H', p.annotation[212] 
    assert p.annotation[232] == 'H', p.annotation[232] 
    assert p.annotation[390:397] == 'EEEEEEE', p.annotation[390:397]
