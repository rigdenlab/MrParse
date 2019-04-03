'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import json
import multiprocessing
import os
import subprocess

from mr_hkl import HklInfo
from mr_search_model import SearchModelFinder
from mrparse.mr_classify import MrClassifier


logger = logging.getLogger(__name__)

HTML_DIR = '/opt/MrParse/html'
HTML_OUT = os.path.join(HTML_DIR, 'mrparse.html')
POLL_TIME = 1


def run(seqin, hklin=None):
    if not (seqin and os.path.isfile(seqin)):
        raise RuntimeError("Cannot find seqin file: %s" % seqin)
    
    logger.info("mr_analyse running with seqin %s", seqin)
     
    # Find homologs and determine properties
    search_model_finder = SearchModelFinder(seqin, hklin=hklin)
    classifier = MrClassifier(seqin=seqin)
    hkl_info = None
    if hklin:
        if not (hklin and os.path.isfile(hklin)):
            raise RuntimeError("Cannot find hklin file: %s" % hklin)
        hkl_info  = HklInfo(hklin)
    
    multip = False
    if multip:
        nproc = 3 if hklin else 2
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
            logger.critical('SearchModelFinder failed: %s' % e)
        try:
            classifier = mrc_result.get()
        except Exception as e:
            logger.critical('MrClassifier failed: %s' % e)
        if hkl_info:
            try:
                hkl_info = hklin_result.get()
            except Exception as e:
                logger.critical('HklInfo failed: %s' % e)
    else:
        search_model_finder()
        classifier()
        if hkl_info:
            hkl_info()

    json_dict = {}
    json_dict.update(search_model_finder.as_dict())
    json_dict.update(classifier.pfam_dict())

    data_dict = {}
    data_dict['pfam'] = json_dict
    if hkl_info:
        data_dict['hkl_info'] = hkl_info.as_dict()
        
    js_data = 'const mrparse_data = %s;\n' % json.dumps(data_dict)
    with open(os.path.join(HTML_DIR, 'mrparse.js'), 'w') as w:
        w.write(js_data)

    # only on Mac OSX
    subprocess.Popen(['open', HTML_OUT])
    return
