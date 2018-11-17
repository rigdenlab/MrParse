'''
Created on 17 Nov 2018

@author: jmht
'''

import set_mrparse_path
import conftest

from mrparse.mr_topcons import transmembrane_prediction


topcons_dir = '/opt/MrParse/data/Q13586/topcons'

predicition = transmembrane_prediction(topcons_dir)


assert len(predicition) == 685
assert predicition[212] == 1.0, predicition[212] 
assert predicition[232] == 1.0, predicition[232] 