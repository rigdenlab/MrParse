#!/usr/bin/env ccp4-python
import set_mrparse_path
import conftest

from mrparse.mr_hkl import HklInfo
from mrparse.mr_sequence import Sequence

def test_2uvo(test_data):
    hkl_info  = HklInfo(test_data.x2uvo_mtz)
    hkl_info()
    assert hkl_info.labels.f == 'FP'
    assert hkl_info.labels.sigf == 'SIGFP'
    assert hkl_info.has_ncs is True
    assert hkl_info.has_twinning is False
    assert hkl_info.has_anisotropy is True


def test_2uvo_seqin(test_data):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    hkl_info = HklInfo(test_data.x2uvo_mtz, seq_info=seq_info)
    assert hkl_info.labels.f == 'FP'
    assert hkl_info.labels.sigf == 'SIGFP'
    assert abs(hkl_info.predicted_solvent_content - 0.508) < 0.01
    assert hkl_info.predicted_ncopies == 4
    assert abs(hkl_info.molecular_weight - 17131) < 0.1

if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])