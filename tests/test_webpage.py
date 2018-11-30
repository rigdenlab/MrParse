'''
Created on 26 Oct 2018

@author: jmht

https://omictools.com/protein-coiled-coil-domain-prediction-category

paircoil2 MIT
deepcoil **


'''
import json
import logging
import os
import pickle
import subprocess
import sys
import set_mrparse_path

from mrparse.mr_sequence import read_fasta
from mrparse.mr_pfam import pfam_region_dict, pfam_dict_from_annotation
from mrparse.mr_classify import MrClassifier
from mrparse.mr_search_model import SearchModelFinder


logging.basicConfig(level=logging.DEBUG)


html_dir = '/opt/MrParse/pfam'
html_out = os.path.join(html_dir, 'mrparse.html')


# seqin = '../data/5u4pA.fasta'
# seqin = '../data/5u4pB.fasta'
# seqin = '../data/O75410.fasta'
seqin = '../data/2uvoA.fasta'
# seqin = '/opt/MrParse/data/Q13586.fasta'

# Find homologs and determine properties
smf = SearchModelFinder(seqin)
regions = smf.find_regions()
smf.find_homologs()
smf.write_html(html_out)

# Get various classification data for graphical display
seqlen = len(read_fasta(seqin))
region_d = pfam_region_dict(regions, seqlen)
 
mrc = MrClassifier()
pfam_data = mrc.get_annotation_pfam_dict(seqin)
 
pfam_data['regions'] = region_d
 
rdata = 'var pfam_json = %s;\n' % json.dumps(pfam_data)
with open(os.path.join(html_dir, 'data.js'), 'w') as w:
    w.write(rdata)

subprocess.Popen(['open', html_out])