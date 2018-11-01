import pickle

import set_mrparse_path
from mrparse.mr_hit import find_hits
from mrparse.mr_region import RegionFinder


def test_region_2uvoA(test_data):
    rfinder = RegionFinder()
    if False:
        seqin = test_data.x2uvoA_fasta
        hits = find_hits(seqin)
        with open('hits.pkl', 'w') as w:
            pickle.dump(hits, w)
    else:
        with open('hits.pkl') as f:
            hits = pickle.load(f)
    regions = rfinder.find_regions_from_hits(hits)
    for r in regions:
        print("GOT REGION %s" % r)

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])