'''
Created on 18 Oct 2018

@author: jmht
'''

import pandas as pd

from mrparse.mr_homolog import get_homologs, calculate_ellg
from mrparse.mr_hit import find_hits
from mr_region import RegionFinder

class SearchModelFinder(object):
    def __init__(self, seqin):
        self.seqin = seqin
        self.hits = None
        self.regions = None
    
    def find_regions(self):
        self.hits = find_hits(self.seqin)
        self.regions = RegionFinder().find_regions_from_hits(self.hits)


    def find_homologs(self, mtz=None):
        self.homologs = get_homologs(self.hits, self.regions)
        if mtz:
            calculate_ellg(self.homologs, mtz)

    def as_dataframe(self):
        homolog_dict = [h.__dict__ for _, h in self.homologs.items()]
        columns = ['name', 'domain', 'range', 'eLLG', 'ncopies', 'molecular_weight', 'rmsd', 'seqid',
                   'frac_scat', 'total_frac_scat', 'total_frac_scat_known', 'pdb']
        df = pd.DataFrame(homolog_dict, columns=columns)
        #df.sort_values('eLLG', inplace=True, ascending=False)
        df.sort_values('domain', inplace=True, ascending=True)
        return df
    
    def as_html(self):
        return self.as_dataframe().to_html(index=False)
