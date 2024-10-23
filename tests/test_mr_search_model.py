#!/usr/bin/env ccp4-python
import set_mrparse_path

import pytest
import logging
from mrparse.mr_search_model import SearchModelFinder
from mrparse.mr_hkl import HklInfo
from mrparse.mr_sequence import Sequence

logging.basicConfig(level=logging.DEBUG)


def test_SearchModelFinder2uvoOnlySequence(test_data):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    smf = SearchModelFinder(seq_info, search_engine="phmmer", pdb_dir=test_data.pdb_dir, phmmer_dblvl="95")
    smf()
    nregions = len(smf.regions)
    assert nregions == 6, f"Incorrect number of regions: {nregions}"
    nhomologs = len(smf.homologs)
    assert nhomologs == 9, f"Incorrect number of homologs: {nhomologs}"
    mw = smf.homologs['1iqb_B'].molecular_weight
    assert abs(mw-8597) < 0.1, f"Incorrect MW: {mw}"


@pytest.mark.skip(reason="Running phaser to calculate eLLG currently takes too long")
def test_SearchModelFinder2uvo(test_data):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hklin = test_data.x2uvo_mtz
    hkl_info = HklInfo(hklin, seq_info=seq_info)
    smf = SearchModelFinder(seq_info, hkl_info=hkl_info, search_engine="phmmer", pdb_dir=test_data.pdb_dir)
    smf()
    nregions = len(smf.regions)
    assert nregions == 6, f"Incorrect number of regions: {nregions}"
    nhomologs = len(smf.homologs)
    assert nhomologs == 9, f"Incorrect number of homologs: {nhomologs}"
    mw = smf.homologs['1iqb_B_1'].molecular_weight
    assert abs(mw - 8597) < 0.1, f"Incorrect MW: {mw}"

if __name__ == '__main__':
    import sys
    pytest.main([__file__] + sys.argv[1:])
