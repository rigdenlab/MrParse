#!/usr/bin/env ccp4-python
import set_mrparse_path
from mrparse.mr_region import RegionFinder


def test_region_2uvoA(get_2uvo_test_hits):
    rfinder = RegionFinder()
    hits = get_2uvo_test_hits
    regions = rfinder.find_regions_from_hits(hits)
    assert len(regions) == 6
    assert regions[2].ID == 3
    assert regions[2].matches[0].name == '1ulm_B_1'

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])