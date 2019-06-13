'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
from mrparse import mr_homolog 
from mrparse import mr_hit
from mr_region import RegionFinder
from mrparse import mr_pfam
from mrparse.mr_util import now

logger = logging.getLogger(__name__)


class SearchModelFinder(object):
    def __init__(self, seq_info, hkl_info=None):
        self.seq_info = seq_info
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
        self.hits = mr_hit.find_hits(self.seq_info)
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

    def pfam_dict(self):
        mr_pfam.add_pfam_dict_to_homologs(self.regions, self.seq_info.nresidues)
        return {'homologs' : [h.json_dict() for h in self.homologs.values()]}
