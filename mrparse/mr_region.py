'''
Created on 18 Oct 2018

@author: jmht

'''
from operator import attrgetter
from mrparse.mr_hit import sort_hits_by_size

class RegionData:
    def __init__(self):
        self.target_name = ""
        self.index = 0
        self.midpoint = 0
        self.extent = 0
        self.matches = []
    
    @property
    def start_stop(self):
        assert self.midpoint >= 0 or self.extent >= 0, "Need non-zero midpoint and extent!"
        half_len = int(self.extent / 2)
        start = self.midpoint - half_len
        stop = self.midpoint + half_len
        return start, stop

    @property
    def id(self):
        return self.index + 1
        
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

    def find_regions_from_hits(self, hits, sort=True):
        """Figure out the regions for the target that have been matched"""
        # Hits need to be sorted from smallest to largest or the domain finding won't work
        if sort:
            hits = sort_hits_by_size(hits, ascending=True)
        target_regions = []
        for hit in hits.values():
#             print "CHECKING HIT %s %s %s" % (hit.name, hit.hit_extent, hit.hit_midpoint)
            self.create_or_update_region(hit, target_regions)
        if sort:
            target_regions = self.sort_regions(target_regions)
        return target_regions
    
    def create_or_update_region(self, hit, target_regions):
        for region in target_regions:
#             print "Checking region %s %s %s" % (region.id, region.extent, region.midpoint)
            if self.hit_within_region(hit, region):
#                 print "WITHIN"
                return self.update_region(hit, region)
        self.add_new_region(hit, target_regions)
        return

    def hit_within_region(self, hit, region, extentTolerance=50, midpointTolerance=20):
#         print "e- %s e+ %s m- %s m+ %s" % (region.extent - extentTolerance,
#                                            region.extent + extentTolerance,
#                                            region.midpoint - midpointTolerance,
#                                            region.midpoint + midpointTolerance)
        if hit.query_extent >= region.extent - extentTolerance and \
            hit.query_extent <= region.extent + extentTolerance and \
            hit.query_midpoint >= region.midpoint - midpointTolerance and \
            hit.query_midpoint <= region.midpoint + midpointTolerance:
            return True
        return False

    def update_region(self, hit, region):
        # Should we update the midpoint and extent of the region?
        region.matches.append(hit)
        hit.region = region
        return
    
    def add_new_region(self, hit, target_regions):
        region = RegionData()
        region.index = len(target_regions)
        region.midpoint = hit.query_midpoint
        region.extent = hit.query_extent
        region.matches.append(hit)
        target_regions.append(region)
        hit.region = region
        return target_regions
    
    def sort_regions(self, regions, ascending=False):
        reverse = not(ascending)
        # Need to think about better ways of sorting - probably store reference to hit in region?
        regions = sorted(regions, key=attrgetter('extent'), reverse=reverse)
        # The matches and ranges also need to be sorted
        for i, r in enumerate(regions):
            r.index = i
            r.matches.reverse()
        return regions
