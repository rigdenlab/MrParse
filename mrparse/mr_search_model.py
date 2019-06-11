'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import pandas as pd

from mrparse import mr_homolog 
from mrparse import mr_hit
from mr_region import RegionFinder
from mrparse import mr_sequence
from mrparse import mr_pfam
from mrparse.mr_util import now

logger = logging.getLogger(__name__)

class SearchModelFinder(object):
    def __init__(self, seqin, hkl_info=None):
        self.seqin = seqin
        self.seqlen = len(mr_sequence.read_fasta(self.seqin))
        self.hkl_info = hkl_info
        self.hits = None
        self.regions = None
        self.homologs = None
    
    def __call__(self):
        """Required so that we can use multiprocessing pool. We need to be able to pickle the object passed
        to the pool and instance methods don't work, so we add the object to the pool and define __call__
        https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/6975654#6975654
        """
        logger.debug('SearchModelFinder started at %s' % now())
        self.find_regions()
        logger.debug('SearchModelFinder regions done at %s' % now())
        self.find_homologs()
        logger.debug('SearchModelFinder homologs done at %s' % now())
        return self
    
    def find_regions(self):
        self.hits = mr_hit.find_hits(self.seqin)
        if not self.hits:
            logger.critical('SearchModelFinder could not find any hits!')
            return None
        self.regions = RegionFinder().find_regions_from_hits(self.hits)
        return self.regions

    def find_homologs(self):
        assert self.hits and self.regions
        self.homologs = mr_homolog.homologs_from_hits(self.hits)
        if self.hkl_info:
            mr_homolog.calculate_ellg(self.homologs, self.hkl_info)
        return self.homologs

    def as_dataframe(self):
        homolog_dict = [h.__dict__ for h in self.homologs.values()]
        columns = ['name', 'domain', 'range', 'eLLG', 'ncopies', 'molecular_weight', 'rmsd', 'seqid',
                   'score', 'frac_scat', 'total_frac_scat', 'total_frac_scat_known', 'pdb_url']
        df = pd.DataFrame(homolog_dict, columns=columns)
        #df.sort_values('eLLG', inplace=True, ascending=False)
        df.sort_values(['domain', 'eLLG'], inplace=True, ascending=True)
        return df
 
    def as_html(self):
        df = self.as_dataframe()
        # Sort number of decimal places
        df = df.round({'molecular_weight': 1,
                       'rmsd' : 2,
                       'score' : 2,
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
        mr_pfam.add_pfam_dict_to_homologs(self.regions, self.seqlen)
        return {'homologs' : [h.json_dict() for h in self.homologs.values()]}
