#!/usr/bin/env ccp4-python
import set_mrparse_path
from mrparse.mr_hit import find_hits, sort_hits_by_size

def Xtest_hit_2uvoA(test_data):
    seqin = test_data.x2uvoA_fasta
    hits = find_hits(seqin)
    hit_names = hits.keys()
    assert hit_names == ['2x3t_C_1', '1ulk_B_1', '1ulm_B_1', '4wp4_A_1', '1iqb_B_1', '5wuz_A_1', '5xdi_A_1', '2lb7_A_1', '1mmc_A_1', '2kus_A_1', '2n1s_A_1'], \
        "Incorrect hit_names: %s" % hit_names

def test_hit_2uvoA_sort_by_size(test_data):
    seqin = test_data.x2uvoA_fasta
    hits = find_hits(seqin)
    hit_names = hits.keys()
    assert hit_names == ['2x3t_C_1', '1ulk_B_1', '1ulm_B_1', '4wp4_A_1', '1iqb_B_1', '5wuz_A_1', '5xdi_A_1', '2lb7_A_1', '1mmc_A_1', '2kus_A_1', '2n1s_A_1'], \
        "Incorrect hit_names: %s" % hit_names
    
    descending = ['2x3t_C_1', '1ulk_B_1', '1ulm_B_1', '1iqb_B_1', '2lb7_A_1', '4wp4_A_1', '5wuz_A_1', '5xdi_A_1', '2kus_A_1', '1mmc_A_1', '2n1s_A_1']
    hits = sort_hits_by_size(hits)
    hit_names = hits.keys()
    assert hit_names == descending, "Incorrect descending: %s" % hit_names
    # 4wp4_A_1 and 5xdi_A_1 are same length
    ascending = [ ['2n1s_A_1', '1mmc_A_1', '2kus_A_1', '5xdi_A_1', '4wp4_A_1', '5wuz_A_1', '2lb7_A_1', '1iqb_B_1', '1ulm_B_1', '1ulk_B_1', '2x3t_C_1'],
                 ['2n1s_A_1', '1mmc_A_1', '2kus_A_1', '5xdi_A_1', '5wuz_A_1', '4wp4_A_1', '2lb7_A_1', '1iqb_B_1', '1ulm_B_1', '1ulk_B_1', '2x3t_C_1']
                ]
    hits = sort_hits_by_size(hits, ascending=True)
    hit_names = hits.keys()
    assert hit_names in ascending, "Incorrect ascending: %s" % hit_names


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])