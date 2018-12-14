'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import json
import multiprocessing
import os
import pickle
import subprocess
import time

from mr_hkl import HklInfo
from mr_search_model import SearchModelFinder
from mrparse.mr_classify import MrClassifier


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

HTML_DIR = '/opt/MrParse/pfam'
POLL_TIME = 1

def write_html(html_out, html_data, template_file='multi_domains_template.html', template_dir='/opt/MrParse/pfam/'):
    from jinja2 import Environment, FileSystemLoader
    env = Environment( 
        loader=FileSystemLoader(template_dir), 
#             autoescape=select_autoescape(['html']) 
        )                                                                                                                                                                          
    template = env.get_template(template_file)
    with open(html_out, 'w') as w:
        w.write(template.render(html_data))
        

# def create_webpage(search_model_finder, classifier, hklinfo=None, html_dir=None):
    


def run(seqin, hklin=None):
    if not os.path.isfile(seqin):
        raise RuntimeError("Cannot find seqin file: %s" % seqin)
     
    # Find homologs and determine properties
    smf = SearchModelFinder(seqin, hklin=hklin)
    mrc = MrClassifier(seqin=seqin)
    if hklin:
        if not os.path.isfile(hklin):
            raise RuntimeError("Cannot find hklin file: %s" % hklin)
        hkl_info  = HklInfo(hklin)
    
    multip = False
    if multip:
        nproc = 3 if hklin else 2
        logger.info("Running on %d processors." % nproc)
        pool = multiprocessing.Pool(nproc)
        smf_result = pool.apply_async(smf)
        mrc_result = pool.apply_async(mrc)
        if hklin:
            hklin_result = pool.apply_async(hkl_info)    
        pool.close()
        logger.debug("Pool waiting")
        pool.join()
        logger.debug("Pool finshed")
        
        # Should raise any child process exceptions
        smf = smf_result.get()
        mrc = mrc_result.get()
        if hklin:
            hkl_info = hklin_result.get()
    else:
        smf()
        mrc()
        if hklin:
            hkl_info()
    
    pfam_data = mrc.pfam_dict()
    pfam_data['regions'] = smf.pfam_dict()
      
    js_data = 'var pfam_json = %s;\n' % json.dumps(pfam_data)
    with open(os.path.join(HTML_DIR, 'data.js'), 'w') as w:
        w.write(js_data)
     
    html_data = {'homolog_table' : smf.as_html()}
    if hklin:
        html_data['hkl_info'] = hkl_info.as_html()
     
    with open(os.path.join(HTML_DIR, 'html_data.pkl'), 'w') as w:
        pickle.dump(html_data, w)

    html_out = os.path.join(HTML_DIR, 'mrparse.html')
    write_html(html_out, html_data)

    # only on Mac OSX
    subprocess.Popen(['open', html_out])
    return