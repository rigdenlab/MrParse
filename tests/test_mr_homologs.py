#!/usr/bin/env ccp4-python
import set_mrparse_path
import pytest
import logging
from mrparse.mr_sequence import Sequence
from mrparse.mr_hkl import HklInfo
from mrparse.mr_homolog import homologs_from_hits, calculate_ellg, ellg_data_from_phaser_log


logging.basicConfig(level=logging.DEBUG)


def test_2uvoA_homologs(test_data, get_2uvo_test_hits):
    hits = get_2uvo_test_hits
    homologs = homologs_from_hits(hits, pdb_dir=test_data.pdb_dir)
    h = homologs['1iqb_B']
    assert h.length == 76
    assert abs(h.molecular_weight - 8597) < 0.001
    assert h.pdb_id == '1iqb'
    assert h.resolution == 0.00
    chain_id = 'B'
    assert h.chain_id == chain_id
    d = h.static_dict
    assert d['chain_id'] == chain_id
    
    
@pytest.mark.skip(reason="Tests using phaser are currently too slow to run")
def test_2uvoA_homologs_ellg(test_data, get_2uvo_test_hits):
    hits = get_2uvo_test_hits
    homologs = homologs_from_hits(hits, pdb_dir=test_data.pdb_dir)
    h = homologs['1iqb_B']
    assert h.length == 77
    assert abs(h.molecular_weight - 8597) < 0.001
    
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hkl_info = HklInfo(test_data.x2uvo_mtz, seq_info=seq_info)
    calculate_ellg(homologs, hkl_info)


def test_phaser_log_ellg(test_data, get_2uvo_test_hits):
    """Test the output from phaser without running it every time"""
    hits = get_2uvo_test_hits
    homologs = homologs_from_hits(hits, pdb_dir=test_data.pdb_dir)
    ellg_data = ellg_data_from_phaser_log(test_data.phaser_log, homologs)
    assert ellg_data['6stq_B'].ellg == 387.3
    assert ellg_data['6stq_B'].length == 166
    assert ellg_data['6stq_B'].molecular_weight == pytest.approx(16194)


if __name__ == '__main__':
    import sys
    pytest.main([__file__] + sys.argv[1:])
