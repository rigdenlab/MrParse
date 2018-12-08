'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import json
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


def run(hklin, seqin):
    
    assert os.path.isfile(seqin)
    if hklin:
        assert os.path.isfile(hklin)
        hkl_info  = HklInfo(hklin)
     
    # Find homologs and determine properties
    smf = SearchModelFinder(seqin)
    smf.find_regions()
    smf.find_homologs(hklin=hklin)
     
    mrc = MrClassifier()
    pfam_data = mrc.get_annotation_pfam_dict(seqin)
    pfam_data['regions'] = smf.get_pfam_dict()
      
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