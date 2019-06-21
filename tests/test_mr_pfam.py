#!/usr/bin/env ccp4-python
'''
Created on 26 Oct 2018

@author: jmht
'''
import set_mrparse_path
from mrparse.mr_sequence import Sequence
from mrparse.mr_homolog import homologs_from_hits
from mrparse.mr_region import RegionFinder
from mrparse.mr_pfam import add_pfam_dict_to_homologs


def test_homologs(test_data, get_2uvo_test_hits):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hits = get_2uvo_test_hits
    regions = RegionFinder().find_regions_from_hits(hits, pdb_download_dir=test_data.pdb_dir)
    homologs = homologs_from_hits(hits)
    add_pfam_dict_to_homologs(regions, seq_info.nresidues)
    
    homolog = homologs.values()[0]
    ali_start = homolog.hit.query_start
    pfam_ali_start = homolog._pfam_json['regions'][0]['aliStart']
    assert ali_start == pfam_ali_start
#     {'regions': [{'startStyle': 'curved', 'end': 79, 'start': 2, 'aliStart': 2, 
#                   'text': '1iqb_B_1', 'endStyle': 'curved', 'colour': '#7f7f3f', 
#                   'aliEnd': 79, 'metadata': 
#                   {'start': 2, 'end': 79, 'description': 'Homolog 1iqb_B_1 from region #2', 
#                    'database': 'PHMMER search'}}], 'length': 171}
    

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
