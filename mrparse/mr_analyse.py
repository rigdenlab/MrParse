'''
Created on 18 Oct 2018

@author: jmht
'''
import logging
import json
import os
import subprocess

from mr_hkl import HklInfo
from mr_search_model import SearchModelFinder
from mrparse.mr_classify import MrClassifier


logging.basicConfig(level=logging.DEBUG)


def write_html(self, html_out, options, template_file='multi_domains_template.html', template_dir='/opt/MrParse/pfam/'):
    from jinja2 import Environment, FileSystemLoader
    env = Environment( 
        loader=FileSystemLoader(template_dir), 
#             autoescape=select_autoescape(['html']) 
        )                                                                                                                                                                          
    template = env.get_template(template_file)
    with open(html_out, 'w') as w:
        w.write(template.render(options))


def run(hklin, seqin):
    html_dir = '/opt/MrParse/pfam'
    html_out = os.path.join(html_dir, 'mrparse.html')
    
    if hklin:
        hkl_info  = HklInfo(hklin)
    
    # Find homologs and determine properties
    smf = SearchModelFinder(seqin)
    smf.find_regions()
    smf.find_homologs(hklin=hklin)
    
    mrc = MrClassifier()
    pfam_data = mrc.get_annotation_pfam_dict(seqin)
    pfam_data['regions'] = smf.get_pfam_dict()
     
    js_data = 'var pfam_json = %s;\n' % json.dumps(pfam_data)
    with open(os.path.join(html_dir, 'data.js'), 'w') as w:
        w.write(js_data)
    
    options = {'homolog_table' : smf.as_html()}
    if hklin:
        options['hkl_info'] = hkl_info.as_html()
    write_html(html_out, options)

    # only on Mac OSX
    subprocess.Popen(['open', html_out])
    return