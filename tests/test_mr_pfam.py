#!/usr/bin/env ccp4-python
"""
Created on 26 Oct 2018

@author: jmht
"""
import set_mrparse_path
from mrparse.mr_sequence import Sequence
from mrparse.mr_homolog import homologs_from_hits
from mrparse.mr_region import RegionFinder
from mrparse.mr_pfam import add_pfam_dict_to_homologs


def test_homologs(test_data, get_2uvo_test_hits):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hits = get_2uvo_test_hits
    RegionFinder().find_regions_from_hits(hits)
    homologs = homologs_from_hits(hits, pdb_dir=test_data.pdb_dir)
    add_pfam_dict_to_homologs(homologs, seq_info.nresidues)

    assert len(homologs) == 14
    homolog = homologs["2x3t_C_1"]
    assert homolog.name == '2x3t_C_1'
    ali_start = homolog.hit.query_start
    pfam_ali_start = homolog._pfam_json['regions'][0]['aliStart']
    assert ali_start == pfam_ali_start


if __name__ == '__main__':
    import sys
    import pytest

    pytest.main([__file__] + sys.argv[1:])
