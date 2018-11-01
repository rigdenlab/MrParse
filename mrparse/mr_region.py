'''
Created on 18 Oct 2018

@author: jmht
'''
from operator import attrgetter
from mrparse.mr_hit import sort_hits_by_size

class RegionData:
    def __init__(self):
        self.targetName = ""
        self.ID = 0
        self.midpoint = 0
        self.extent = 0
        self.matches = []
        self.ranges = []
    
    @property
    def start_stop(self):
        assert self.midpoint >= 0 or self.extent >= 0, "Need non-zero midpoint and extent!"
        half_len = int(self.extent / 2)
        start = self.midpoint - half_len
        stop = self.midpoint + half_len
        return start, stop
        
    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str


class RegionFinder(object):
    def __init__(self):
        pass

    def find_regions_from_hits(self, hits):
        """Figure out the regions for the target that have been matched"""
        # Hits need to be sorted from smallest to largest or the domain finding won't work
        hits = sort_hits_by_size(hits, ascending=True)
        targetRegions = []
        for hit in hits.values():
            self.create_or_update_region(hit, targetRegions)
        return self.sort_regions(targetRegions)
    
    def create_or_update_region(self, hit, targetRegions):
        for region in targetRegions:
            if self.hit_within_region(hit, region):
                return self.update_region(hit, region)
        self.add_new_region(hit, targetRegions)
        return

    def hit_within_region(self, hit, region, extentTolerance=50, midpointTolerance=20):
        if hit.tarExtent >= region.extent - extentTolerance and \
            hit.tarExtent <= region.extent + extentTolerance and \
            hit.tarMidpoint >= region.midpoint - midpointTolerance and \
            hit.tarMidpoint <= region.midpoint + midpointTolerance:
            return True
        return False

    def update_region(self, hit, region):
        # Should we update the midpoint and extent of the region?
        region.matches.append(hit.name)
        region.ranges.append(hit.tarRange)
        return
    
    def add_new_region(self, hit, targetRegions):
        region = RegionData()
        region.ID = len(targetRegions) + 1
        region.midpoint = hit.tarMidpoint
        region.extent = hit.tarExtent
        region.matches.append(hit.name)
        region.ranges.append(hit.tarRange)
        targetRegions.append(region)
        return targetRegions
    
    def sort_regions(self, regions, ascending=False):
        reverse = not(ascending)
        # Need to think about better ways of sorting - probably store reference to hit in region?
        regions = sorted(regions, key=attrgetter('extent'), reverse=reverse)
        # The matches and ranges also need to be sorted
        for r in regions:
            r.matches.reverse()
            r.ranges.reverse()
        return regions
