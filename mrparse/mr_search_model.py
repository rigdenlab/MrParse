'''
Created on 18 Oct 2018

@author: jmht
'''
import json
import pandas as pd

from mrparse.mr_homolog import get_homologs, calculate_ellg
from mrparse.mr_hit import find_hits
from mr_region import RegionFinder
from mrparse.mr_sequence import read_fasta
from mrparse.mr_pfam import pfam_region_dict

class SearchModelFinder(object):
    def __init__(self, seqin, hklin=None):
        self.seqin = seqin
        self.hklin = hklin
        self.hits = None
        self.regions = None
    
    def execute(self, queue=None):
        self.find_regions()
        self.find_homologs()
        if queue:
            queue.put(self)
    
    def find_regions(self):
        self.hits = find_hits(self.seqin)
        self.regions = RegionFinder().find_regions_from_hits(self.hits)
        return self.regions

    def find_homologs(self):
        self.homologs = get_homologs(self.hits, self.regions)
        if self.hklin:
            calculate_ellg(self.homologs, self.hklin)
        return self.homologs

    def as_dataframe(self):
        homolog_dict = [h.__dict__ for _, h in self.homologs.items()]
        columns = ['name', 'domain', 'range', 'eLLG', 'ncopies', 'molecular_weight', 'rmsd', 'seqid',
                   'frac_scat', 'total_frac_scat', 'total_frac_scat_known', 'pdb_url']
        df = pd.DataFrame(homolog_dict, columns=columns)
        #df.sort_values('eLLG', inplace=True, ascending=False)
        df.sort_values('domain', inplace=True, ascending=True)
        return df
    
    def as_html(self):
        df = self.as_dataframe()
        # Sort number of decimal places
        df = df.round({'molecular_weight': 1,
                       'rmsd' : 2,
                       'frac_scat' : 2,
                       'total_frac_scat' : 2,
                       'total_frac_scat_known' : 2})
        # Convert pdb_urls to html anchors
        def name_and_url_to_anchor(x):
            name = x[0]
            pdb = name.split('_')[0]
            url = x[1]
            return  "<a href=\"%s\">%s</a>" % (url, pdb)
        df['pdb_url'] = df[['name','pdb_url']].apply(name_and_url_to_anchor, axis=1)
        pd.set_option('display.max_colwidth', -1)
        return df.to_html(index=False, escape=False)

    def pfam_dict(self):
        assert self.regions and self.seqin
        seqlen = len(read_fasta(self.seqin))
        return pfam_region_dict(self.regions, seqlen)