#!/usr/bin/env ccp4-python
import set_mrparse_path
from mrparse.mr_region import RegionFinder


def test_region_2uvoA(get_2uvo_test_hits):
    rfinder = RegionFinder()
    hits = get_2uvo_test_hits
    regions = rfinder.find_regions_from_hits(hits)
#     for h in hits.values():
#         #print("GOT HIT {} : {} : {} : {}".format(h.name, h.tar_range, h.tar_extent, h.tar_midpoint))
#         print("GOT HIT {} : {} : {} : {}".format(h.name, h.query_range, h.query_extent, h.query_midpoint))
#     for i, r in enumerate(regions):
#         print("{} matches: {}".format(i, [m.name for m in r.matches]))
    assert len(regions) == 5
    assert regions[2].id == 3
    assert regions[2].matches[0].name == '1iqb_B_1'

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])