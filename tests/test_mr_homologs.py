#!/usr/bin/env ccp4-python
import set_mrparse_path
import pytest
from mrparse.mr_sequence import Sequence
from mrparse.mr_hkl import HklInfo
from mrparse.mr_homolog import homologs_from_hits, calculate_ellg, ellg_data_from_phaser_log


def test_2uvoA_homologs(get_2uvo_test_hits):
    hits = get_2uvo_test_hits
    homologs = homologs_from_hits(hits)
    h = homologs['1iqb_B_1']
    assert h.length == 77
    assert abs(h.molecular_weight - 8750.78) < 0.001
    

@pytest.mark.skip(reason="Tests using phaser are currently too slow to run")
def test_2uvoA_homologs_ellg(test_data, get_2uvo_test_hits):
    hits = get_2uvo_test_hits
    homologs = homologs_from_hits(hits)
    h = homologs['1iqb_B_1']
    assert h.length == 77
    assert abs(h.molecular_weight - 8750.78) < 0.001
    
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hkl_info = HklInfo(test_data.x2uvo_mtz, seq_info=seq_info)
    calculate_ellg(homologs, hkl_info)


def test_phaser_log_ellg(test_data, get_2uvo_test_hits):
    """Test the output from phaser without running it every time"""
    hits = get_2uvo_test_hits
    homologs = homologs_from_hits(hits)
    ellg_data = ellg_data_from_phaser_log(test_data.phaser_log, homologs)
    assert ellg_data['2x3t_C_1'].ellg == 1378.4
    assert ellg_data['2x3t_C_1'].length == 168
    assert ellg_data['2x3t_C_1'].molecular_weight == pytest.approx(16678.6458)

if __name__ == '__main__':
    import sys
    pytest.main([__file__] + sys.argv[1:])