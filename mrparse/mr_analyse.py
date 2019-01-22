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

from mr_hkl import HklInfo
from mr_search_model import SearchModelFinder
from mrparse.mr_classify import MrClassifier


logger = logging.getLogger(__name__)

HTML_DIR = '/opt/MrParse/html'
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


    logger.critical("!!! HACK FOR JSON!!!!")
    search_model_finder.as_json()
    with open(os.path.join(HTML_DIR, 'homologs.js'), 'w') as w:
        w.write(search_model_finder.as_json())
    
    html_out = create_webpage(search_model_finder,
                              classifier,
                              hkl_info=hkl_info,
                              html_dir=HTML_DIR)
    # only on Mac OSX
    subprocess.Popen(['open', html_out])
    return


def create_webpage(search_model_finder,
                   classifier,
                   hkl_info=None,
                   html_dir=None,
                   html_filename='mrparse.html'):
    
    pfam_data = classifier.pfam_dict()
    pfam_data['regions'] = search_model_finder.pfam_dict()
      
    js_data = 'var pfam_json = %s;\n' % json.dumps(pfam_data)
    with open(os.path.join(html_dir, 'data.js'), 'w') as w:
        w.write(js_data)
     
    html_data = {'homolog_table' : search_model_finder.as_html()}
    if hkl_info:
        html_data['hkl_info'] = hkl_info.as_html()
     
    with open(os.path.join(html_dir, 'html_data.pkl'), 'w') as w:
        pickle.dump(html_data, w)

    html_out = os.path.join(html_dir, html_filename)
    write_html(html_out, html_data)
    return html_out

def write_html(html_out, html_data, template_file='multi_domains_template.html', template_dir='/opt/MrParse/pfam/'):
    from jinja2 import Environment, FileSystemLoader
    env = Environment( 
        loader=FileSystemLoader(template_dir), 
#             autoescape=select_autoescape(['html']) 
        )                                                                                                                                                                          
    template = env.get_template(template_file)
    with open(html_out, 'w') as w:
        w.write(template.render(html_data))
