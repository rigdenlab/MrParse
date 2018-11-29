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
import sys
import set_mrparse_path

from mrparse.mr_sequence import read_fasta
from mrparse.mr_pfam import pfam_region_dict, pfam_classification_dict
from mrparse.mr_classify import get_annotation
from mrparse.mr_annotation import get_annotation_chunks
from mrparse.mr_jpred import JPred
from mrparse.mr_search_model import SearchModelFinder


logging.basicConfig(level=logging.DEBUG)


html_dir = '/opt/MrParse/pfam'
# seqin = '../data/5u4pA.fasta'
# seqin = '../data/5u4pB.fasta'
# seqin = '../data/O75410.fasta'
seqin = '../data/2uvoA.fasta'
# seqin = '../O75410_800/O75410_800.fasta'
# seqin = '/opt/MrParse/data/Q13586.fasta'

find_homologs = True
if find_homologs:
    smf = SearchModelFinder(seqin)
    smf.find_regions()
    smf.find_homologs()
    regions = smf.regions
    html = json.dumps(smf.as_html())
    with open(os.path.join(html_dir, 'homologs.html.js'), 'w ') as w:
        s = """
        function homologs()
{{
        document.write({});
}}
""".format(html)
        w.write(s)
else:
    from mrparse.mr_region import RegionFinder
    from mrparse.mr_hit import find_hits
    rfinder = RegionFinder()
    hits = find_hits(seqin)
#     with open('hits.pkl', 'w') as w:
#         pickle.dump(hits, w)
    regions = rfinder.find_regions_from_hits(hits)
    with open('regions.pkl', 'w') as w:
        pickle.dump(regions, w)
# else:
#     with open('regions.pkl') as regions_fh:
#         regions = pickle.load(regions_fh)
 
seqlen = len(read_fasta(seqin))
region_d = pfam_region_dict(regions, seqlen)

annotation = get_annotation(seqin)
annotation_chunks = get_annotation_chunks(annotation)
class_d = pfam_classification_dict(annotation_chunks, len(annotation))
 
sspred = JPred().secondary_structure_prediction(seqin=seqin)
sspred_chunks = get_annotation_chunks(sspred)
sspred_d = pfam_classification_dict(sspred_chunks, len(sspred))

pfam_data = {}
pfam_data['ss_pred'] = sspred_d
pfam_data['classification'] = class_d
pfam_data['regions'] = region_d

rdata = 'console.log("LOADED DATA");\nvar pfam_json = %s;\n' % json.dumps(pfam_data)
with open(os.path.join(html_dir, 'data.js'), 'w') as w:
    w.write(rdata)
