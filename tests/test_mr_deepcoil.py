#!/usr/bin/env ccp4-python
import set_mrparse_path
import pytest
import configparser
import logging
import os
from mrparse.mr_sequence import Sequence
from mrparse.mr_deepcoil import CCPred

logging.basicConfig(level=logging.DEBUG)

config_file = os.path.join(os.environ["CCP4"], "share", "mrparse", "data", "mrparse.config")
config = configparser.ConfigParser()
config.read(config_file)
exe_dict = dict(config.items("Executables"))


@pytest.mark.skip(reason="Deepcoil exe needs to be available.")
def test_run(test_data):
    seq_info = Sequence(test_data.x2uvoA_fasta)
    cc = CCPred(seq_info, exe_dict['deepcoil_exe'])
    cc.get_prediction()
    annotation = cc.prediction
    assert len(annotation) == 171
    assert annotation.scores[0] > 0.01 # first few are always non-zero


if __name__ == '__main__':
    import sys
    import pytest
    pytest.main([__file__] + sys.argv[1:])
