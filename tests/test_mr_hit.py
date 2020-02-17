#!/usr/bin/env ccp4-python
import set_mrparse_path
from mrparse.mr_sequence import Sequence
from mrparse.mr_hit import find_hits, sort_hits_by_size


def test_hit_2uvoA(test_data):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hits = find_hits(seq_info)
    hit_names = hits.keys()
    ref_hits = ['2uvo_F_1']
    # hits could change as new sequences are added so we just check for inclusion
    assert set(ref_hits).issubset(set(hit_names)), "Incorrect hit_names: %s" % hit_names
        

def test_hit_2uvoA_sort_by_size(test_data, get_2uvo_test_hits):
    """Hit tar_extent:
    '2x3t_C_1', 166
    '1ulk_B_1', 117
    '1eis_A_1', 82
    '1iqb_B_1', 82
    '1ulm_B_1', 79
    '2lb7_A_1', 41
    '4wp4_A_1', 39
    '5wuz_A_1', 39
    '4mpi_A_1', 37
    '5xdi_A_1', 36
    '2kus_A_1', 29
    '1mmc_A_1', 23
    '1zuv_A_1', 22
    '2n1s_A_1', 22
    """
    
    # Use hits with know file so we always get the same results
    hits = get_2uvo_test_hits
    hit_names = hits.keys()
    assert hit_names == ['2x3t_C_1', '1ulk_B_1', '1ulm_B_1', '4wp4_A_1', '1eis_A_1', '1iqb_B_1', '4mpi_A_1',
                         '5wuz_A_1', '5xdi_A_1', '2lb7_A_1', '1zuv_A_1', '1mmc_A_1', '2kus_A_1', '2n1s_A_1'], \
                         "Incorrect hit_names: %s" % hit_names
    # Some have same sizes so can only check some
    hits = sort_hits_by_size(hits)
    hit_names = hits.keys()
    name = '2x3t_C_1'
    assert hit_names.index(name) == 0, "Incorrect descending for: %s" % name
    name = '2lb7_A_1'
    assert hit_names.index(name) == 5, "Incorrect descending for: %s" % name
    name = '1mmc_A_1'
    assert hit_names.index(name) == 11, "Incorrect descending for: %s" % name

    hits = sort_hits_by_size(hits, ascending=True)
    hit_names = hits.keys()
    name = '2x3t_C_1'
    assert hit_names.index(name) == 13, "Incorrect ascending for: %s" % name
    name = '2lb7_A_1'
    assert hit_names.index(name) == 8, "Incorrect ascending for: %s" % name
    name = '1mmc_A_1'
    assert hit_names.index(name) == 2, "Incorrect ascending for: %s" % name


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])