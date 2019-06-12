#!/usr/bin/env ccp4-python
import set_mrparse_path
from mrparse.mr_homolog import homologs_from_hits


def test_region_2uvoA(test_data, get_2uvo_test_hits):
    hits = get_2uvo_test_hits
    homologs = homologs_from_hits(hits)
    h = homologs['1iqb_B_1']
    assert h.length == 77
    assert abs(h.molecular_weight - 8750.78) < 0.001

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])