import os
import pytest
from collections import namedtuple

from mrparse.mr_hit import _find_hits
import data_constants


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
            'x2uvo_mtz' : os.path.join(data_dir, '2uvo_pdbredo.mtz')
        }
    nt = namedtuple('TestData', d.keys())
    return nt(**d)

@pytest.fixture(scope="session")
def get_2uvo_test_hits():
    phmmer_logfile = 'phmmer_test.log'
    with open(phmmer_logfile, 'w') as w:
        w.write(data_constants.PHMMER_LOG_TXT)
    hits = _find_hits(logfile=phmmer_logfile,
                      searchio_type='hmmer3-text',
                      target_sequence=data_constants.TWOUVO_SEQ)
    os.unlink(phmmer_logfile)
    return hits
     
    
    
