'''
Created on 18 Oct 2018

@author: jmht
'''

import pandas as pd

from mrparse.mr_homolog import get_homologs, ellg_data_from_phaser_log, calculate_ellg
from mrparse.mr_hit import find_hits
from mr_region import RegionFinder

class SearchModelFinder(object):
    def __init__(self, seqin):
        self.seqin = seqin
        self.find_homologs()
    
    def find_homologs(self, mock=True):
#         import pickle
#         with open('../tests/homologs.pkl') as f:
#             self.homologs = pickle.load(f)
#             return
        hits = find_hits(self.seqin)
        domains = RegionFinder().find_regions_from_hits(hits)
        homologs = get_homologs(hits, domains)
        if mock:
            ellg_data_from_phaser_log('phaser1.log', homologs)
        else:
            calculate_ellg(homologs)
        self.homologs = homologs
        
    def as_dataframe(self):
        homolog_dict = [h.__dict__ for _, h in self.homologs.items()]
        columns = ['name', 'domain', 'range', 'eLLG', 'ncopies', 'molecular_weight', 'rmsd', 'seqid',
                   'frac_scat', 'total_frac_scat', 'total_frac_scat_known', 'pdb']
        df = pd.DataFrame(homolog_dict, columns=columns)
        df.sort_values('eLLG', inplace=True, ascending=False)
        return df
    
    def as_html(self):
        return self.as_dataframe().to_html(index=False)
