import set_mrparse_path
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
            'x2uvoA_fasta': os.path.join(data_dir, '2uvoA.fasta'),
            'x2uvo_mtz': os.path.join(data_dir, '2uvo_pdbredo.mtz'),
            'Q13586_fasta': os.path.join(data_dir, 'Q13586.fasta'),
            'topcons_output': os.path.join(data_dir, 'query.result.txt'),
            'jpred_concise': os.path.join(data_dir, 'jp_dSDGWMs.concise'),
            'jpred_tgz': os.path.join(data_dir, 'jp_5or_zBY.tar.gz'),
            'phaser_log': os.path.join(data_dir, 'phaser1.log'),
            'x5hxg_fasta': os.path.join(data_dir, '5hxg.fasta'),
            'pdb_dir': os.path.join(data_dir, 'pdbs')
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


@pytest.fixture(scope="session")
def get_2uvo_alphafold_test_hits():
    phmmer_logfile = 'phmmer_test.log'
    with open(phmmer_logfile, 'w') as w:
        w.write(data_constants.PHMMER_AF_LOG_TXT)
    hits = _find_hits(logfile=phmmer_logfile,
                      searchio_type='hmmer3-text',
                      target_sequence=data_constants.TWOUVO_SEQ)
    os.unlink(phmmer_logfile)
    return hits
