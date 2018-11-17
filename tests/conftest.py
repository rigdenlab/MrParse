import os
import pytest
from collections import namedtuple

# if not "CCP4" in sorted(os.environ.keys()):
#     raise RuntimeError('CCP4 not found')
# CCP4 = os.environ["CCP4"]

# Pytest fixtures
@pytest.fixture(scope="session")
def test_data():
    """Return a namedtuple object with the paths to all the data files we require.

    Returns
    -------
    :obj:`namedtuple <collections.namedtuple>`
        object with paths to required data files
    """
    data_dir = '../data'
    d = {
        'x2uvoA_fasta' : os.path.join(data_dir, '2uvoA.fasta'),
        'x2uvo_mtz' : os.path.join(data_dir, '2uvo_pdbredo.mtz'),
#         'toxd_fasta' : os.path.join(CCP4, "examples", "toxd", "toxd.seq"),
#         'toxd_pdb' : os.path.join(CCP4, "examples", "toxd", "toxd.pdb"),
#         'toxd_mtz' : os.path.join(CCP4, "examples", "toxd", "toxd.mtz"),
#         'x3a22_pdb' : os.path.join(CCP4, "examples", "data", "3a22.pdb"),
#         'rnase_fasta' : os.path.join(CCP4, "examples", "rnase", "rnase.seq"),
#         'rnase_mtz' : os.path.join(CCP4, "examples", "rnase", "rnase18.mtz"),
#         'hypf_mtz' : os.path.join(CCP4, "examples", "mr_tutorial_2006", "data", "hypF", "hypF-1gxu-1gxt-HG_scaleit1.mtz"),
        }
    nt = namedtuple('TestData', d.keys())
    return nt(**d)
