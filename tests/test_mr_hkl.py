'''
Created on 17 Nov 2018

@author: jmht
'''
import set_mrparse_path
import conftest

from mrparse.mr_hkl import HklInfo

def test_2uvo(test_data):
    hkl_info  = HklInfo(test_data.x2uvo_mtz)
    assert hkl_info.labels.f == 'FP'
#     print(hkl_info.labels.fplus)
#     print(hkl_info.labels.fminus)
#     print(hkl_info.labels.dano)
    assert hkl_info.has_ncs is True
    assert hkl_info.has_twinning is False
    assert hkl_info.has_anisotropy is True
    
    
def test_foo():
    hklin = '/opt/MrParse/svetlana/AxMncat_P21.mtz'
    hkl_info  = HklInfo(hklin)
    print(hkl_info.labels.f)
    print(hkl_info.labels.fplus)
    print(hkl_info.labels.fminus)
    print(hkl_info.labels.dano)
    print(hkl_info.has_ncs)
    print(hkl_info.has_twinning)
    print(hkl_info.has_anisotropy)
