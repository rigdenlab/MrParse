import pickle

import set_mrparse_path
from mrparse.mr_hit import find_hits
from mrparse.mr_region import RegionFinder


def test_region_2uvoA(test_data):
    rfinder = RegionFinder()
    seqin = test_data.x2uvoA_fasta
    hits = find_hits(seqin)
    regions = rfinder.find_regions_from_hits(hits)
    assert len(regions) == 5
    assert regions[2].ID == 3
    assert regions[2].matches == ['5wuz_A_1', '4wp4_A_1']

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])