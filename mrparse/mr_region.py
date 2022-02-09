"""
Created on 18 Oct 2018

@author: jmht
"""

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
        indent = "  "
        out_str = f"Class: {self.__class__}\nData:\n"
        for a in sorted(attrs):
            out_str += indent + f"{a} : {self.__dict__[a]}\n"
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
            self.create_or_update_region(hit, target_regions)
        if sort:
            target_regions = self.sort_regions(target_regions)
        return target_regions

    def create_or_update_region(self, hit, target_regions):
        for region in target_regions:
            if self.hit_within_region(hit, region):
                return self.update_region(hit, region)
        self.add_new_region(hit, target_regions)
        return

    @staticmethod
    def hit_within_region(hit, region, extent_tolerance=50, midpoint_tolerance=20):
        if hit.query_extent >= region.extent - extent_tolerance:
            if hit.query_extent <= region.extent + extent_tolerance:
                if hit.query_midpoint >= region.midpoint - midpoint_tolerance:
                    if hit.query_midpoint <= region.midpoint + midpoint_tolerance:
                        return True
        return False

    @staticmethod
    def update_region(hit, region):
        # Should we update the midpoint and extent of the region?
        region.matches.append(hit)
        hit.region = region
        return

    @staticmethod
    def add_new_region(hit, target_regions):
        region = RegionData()
        region.index = len(target_regions)
        region.midpoint = hit.query_midpoint
        region.extent = hit.query_extent
        region.matches.append(hit)
        target_regions.append(region)
        hit.region = region
        return target_regions

    @staticmethod
    def sort_regions(regions, ascending=False):
        reverse = not ascending
        # Need to think about better ways of sorting - probably store reference to hit in region?
        regions = sorted(regions, key=attrgetter('extent'), reverse=reverse)

        # The matches and ranges also need to be sorted
        for i, region in enumerate(regions):
            region.index = i
            region.matches.reverse()
        return regions
