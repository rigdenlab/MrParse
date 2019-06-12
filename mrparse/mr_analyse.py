'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import json
import multiprocessing
import os
import subprocess
import sys

from mr_hkl import HklInfo
from mr_search_model import SearchModelFinder
from mr_sequence import Sequence
from mr_classify import MrClassifier


logger = logging.getLogger(__name__)

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
HTML_DIR= os.path.join(THIS_DIR, '../html')
HTML_OUT = os.path.join(HTML_DIR, 'mrparse.html')
POLL_TIME = 1


def run(seqin, hklin=None, run_serial=False):
    if not (seqin and os.path.isfile(seqin)):
        raise RuntimeError("Cannot find seqin file: %s" % seqin)
    logger.info("mr_analyse running with seqin %s", os.path.abspath(seqin))
    seq_info = Sequence(seqin)
    hkl_info = None
    if hklin:
        if not os.path.isfile(hklin):
            raise RuntimeError("Cannot find hklin file: %s" % hklin)
        logger.info("mr_analyse running with hklin %s", os.path.abspath(hklin))
        hkl_info  = HklInfo(hklin, seq_info=seq_info)
    search_model_finder = SearchModelFinder(seq_info, hkl_info=hkl_info)
    classifier = MrClassifier(seq_info=seq_info)
    if run_serial:
        try:
            search_model_finder()
        except Exception as e:
            logger.exception('SearchModelFinder failed: %s' % e)
        try:
            classifier()
        except Exception as e:
            logger.exception('MrClassifier failed: %s' % e)
        if hkl_info:
            try:
                hkl_info()
            except Exception as e:
                logger.exception('HklInfo failed: %s' % e)
    else:
        nproc = 3 if hkl_info else 2
        logger.info("Running on %d processors." % nproc)
        pool = multiprocessing.Pool(nproc)
        smf_result = pool.apply_async(search_model_finder)
        mrc_result = pool.apply_async(classifier)
        if hkl_info:
            hklin_result = pool.apply_async(hkl_info)    
        pool.close()
        logger.debug("Pool waiting")
        pool.join()
        logger.debug("Pool finished")
        try:
            search_model_finder = smf_result.get()
        except Exception as e:
            logger.exception('SearchModelFinder failed: %s' % e)
        try:
            classifier = mrc_result.get()
        except Exception as e:
            logger.exception('MrClassifier failed: %s' % e)
        if hkl_info:
            try:
                hkl_info = hklin_result.get()
            except Exception as e:
                logger.exception('HklInfo failed: %s' % e)

    pfam_dict = {}
    pfam_dict.update(search_model_finder.pfam_dict())
    pfam_dict.update(classifier.pfam_dict())
    data_dict = {}
    data_dict['pfam'] = pfam_dict
    if hkl_info:
        data_dict['hkl_info'] = hkl_info.as_dict()
        
    js_data = 'const mrparse_data = %s;\n' % json.dumps(data_dict)
    with open(os.path.join(HTML_DIR, 'mrparse.js'), 'w') as w:
        w.write(js_data)

    # Display results in browser
    opencmd = None
    if sys.platform.lower().startswith('linux'):
        opencmd = 'xdg-open'
    elif sys.platform.lower().startswith('darwin'):
        opencmd = 'open'
    if opencmd:
        subprocess.Popen([opencmd, HTML_OUT])
    return 0
