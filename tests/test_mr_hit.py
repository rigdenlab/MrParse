#!/usr/bin/env ccp4-python
import set_mrparse_path
from mrparse.mr_sequence import Sequence
from mrparse.mr_hit import find_hits, sort_hits_by_size


def test_hit_2uvoA(test_data):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hits = find_hits(seq_info)
    hit_names = hits.keys()
    ref_hits = ['2uvo_F']
    # hits could change as new sequences are added so we just check for inclusion
    assert set(ref_hits).issubset(set(hit_names)), f"Incorrect hit_names: {hit_names}"
        

def test_hit_2uvoA_sort_by_size(test_data, get_2uvo_test_hits):
    """Hit tar_extent:
    2uvo_F: 2-171
    6stq_B: 1-170
    1ulk_B: 201-326
    1uha_A: 1-82
    4wp4_A: 1-43
    1en2_A: 2-86
    1iqb_B: 2-89
    4mpi_B: -1-43
    2dkv_A: 33-330
    """
    
    # Use hits with know file so we always get the same results
    hits = get_2uvo_test_hits
    hit_names = [*hits]

    assert hit_names == ['2uvo_F', '6stq_B', '1ulk_B', '1uha_A', '4wp4_A', '1en2_A', '1iqb_B', '4mpi_B', '2dkv_A'], \
                         f"Incorrect hit_names: {hit_names}"
    # Some have same sizes so can only check some
    hits = sort_hits_by_size(hits)
    hit_names = [*hits]

    name = '2uvo_F'
    assert hit_names.index(name) == 0, f"Incorrect descending for: {name}"
    name = '1en2_A'
    assert hit_names.index(name) == 4, f"Incorrect descending for: {name}"
    name = '2dkv_A'
    assert hit_names.index(name) == 8, f"Incorrect descending for: {name}"

    hits = sort_hits_by_size(hits, ascending=True)
    hit_names = [*hits]

    name = '2uvo_F'
    assert hit_names.index(name) == 8, f"Incorrect ascending for: {name}"
    name = '1ulk_B'
    assert hit_names.index(name) == 6, f"Incorrect ascending for: {name}"
    name = '4mpi_B'
    assert hit_names.index(name) == 2, f"Incorrect ascending for: {name}"


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])