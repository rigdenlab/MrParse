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


logging.basicConfig(level=logging.DEBUG)

HTML_DIR = '/opt/MrParse/pfam'


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
    if hklin:
        if not os.path.isfile(hklin):
            raise RuntimeError("Cannot find hklin file: %s" % hklin)
        hkl_info  = HklInfo(hklin)
     
    # Find homologs and determine properties
    smf = SearchModelFinder(seqin, hklin=hklin)
    mrc = MrClassifier(seqin=seqin)
    
    multip = True
    if multip:
        queue = multiprocessing.Queue()
        p1 = multiprocessing.Process(target=smf.execute, args=(queue,))
        p2 = multiprocessing.Process(target=mrc.execute, args=(queue,))
        p1.start()
        p2.start()
        if hklin:
            p3 = multiprocessing.Process(target=hkl_info.execute, args=(queue,))
            p3.start()
        p1.join()
        p2.join()
        if hklin:
            p3.join()
        while not queue.empty():
            obj = queue.get()
            if isinstance(obj, SearchModelFinder):
                smf = obj
            elif isinstance(obj, MrClassifier):
                mrc = obj
            elif isinstance(obj, HklInfo):
                hkl_info = obj
    else:
        smf.execute()
        mrc.execute()
        if hklin:
            hkl_info.execute()
    
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