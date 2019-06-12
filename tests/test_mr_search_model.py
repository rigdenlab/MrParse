#!/usr/bin/env ccp4-python
import set_mrparse_path

import logging
from mrparse.mr_search_model import SearchModelFinder
from mrparse.mr_hkl import HklInfo


logging.basicConfig(level=logging.DEBUG)


def Xtest_SearchModelFinder2uvoOnlySequence(test_data):
    seqin = test_data.x2uvoA_fasta
    smf = SearchModelFinder(seqin)
    smf()
    nregions = len(smf.regions)
    assert nregions == 6, "Incorrect number of regions: {}".format(nregions)
    nhomologs = len(smf.homologs)
    assert nhomologs == 14, "Incorrect number of homologs: {}".format(nhomologs)
    assert abs(smf.homologs['1iqb_B_1'].molecular_weight-8750.8) < 0.1, \
    "Incorrect MW: {}".format(smf.homologs['1iqb_B_1'].molecular_weight)
#     print(smf.as_html())
    
def test_SearchModelFinder2uvo(test_data):
    seqin = test_data.x2uvoA_fasta
    hklin = test_data.x2uvo_mtz
    hkl_info = HklInfo(hklin)
    smf = SearchModelFinder(seqin, hkl_info=hkl_info)
    smf()
    nregions = len(smf.regions)
    assert nregions == 6, "Incorrect number of regions: {}".format(nregions)
    nhomologs = len(smf.homologs)
    assert nhomologs == 14, "Incorrect number of homologs: {}".format(nhomologs)
    assert abs(smf.homologs['1iqb_B_1'].molecular_weight-8750.8) < 0.1, \
    "Incorrect MW: {}".format(smf.homologs['1iqb_B_1'].molecular_weight)



if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])

# rfinder = RegionFinder()
# seqin = '../data/2uvoA.fasta'
# hits = find_hits(seqin)
# with open('hits.pkl', 'w') as w:
#     pickle.dump(hits, w)
# regions = rfinder.find_regions_from_hits(hits)
# with open('regions.pkl', 'w') as w:
#     pickle.dump(regions, w)
    
# homologs = get_homologs(hits, regions)
# with open('homologs.pkl', 'w') as w:
#     pickle.dump(homologs, w)
 
# with open('homologs1.pkl') as f:
#     homologs = pickle.load(f) 
# ellg_data_from_phaser_log('phaser1.log', homologs)
 
# calculate_ellg(homologs)
# with open('homologs.pkl', 'w') as w:
#     pickle.dump(homologs, w)
        
# with open('homologs.pkl') as f:
#     homologs = pickle.load(f)
# for h, d in homologs.items():
#     print h, d

# ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py
