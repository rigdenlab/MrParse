'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest

from mrparse.mr_jpred import secondary_structure_prediction


jpred_dir = '/opt/MrParse/data/Q13586/jpred'

predicition = secondary_structure_prediction(jpred_dir)


assert len(predicition) == 685
assert predicition[212] == 'H', predicition[212] 
assert predicition[232] == 'H', predicition[232] 
