'''
Created on 26 Oct 2018

@author: jmht

https://omictools.com/protein-coiled-coil-domain-prediction-category

paircoil2 MIT
deepcoil **


'''
import json
import os
import pickle
import sys
import set_mrparse_path
from mrparse.mr_sequence import read_fasta
from mrparse.mr_classify import MrClassifier
from mrparse.mr_pfam import pfam_region_dict

from mrparse.mr_search_model import SearchModelFinder


html_dir = '/opt/MrParse/pfam'
# seqin = '../data/5u4pA.fasta'
# seqin = '../data/5u4pB.fasta'
# seqin = '../data/O75410.fasta'
# seqin = '../data/2uvoA.fasta'
# seqin = '../data/O75410.fasta'
# seqin = '../O75410_800/O75410_800.fasta'
seqin = '/opt/MrParse/data/Q13586.fasta'

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
        document.write("{}");
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
region_data = pfam_region_dict(regions, seqlen)

# class_data = coiled_coil_pfam_dict(seqin)
# start = 50
# stop = 500
# disordered_region = { 'startStyle': "straight",
#               'endStyle': "straight",
#               'start': start,
#               'end': stop,
#               'aliStart': start,
#               'aliEnd': stop,
#               'colour': "#d3d3d3",
#               'text': 'Disordered',
#               'metadata' : { "description" : "Disordered region #%d" % 1,
#                              "database" : "From Jens",
#                              "start" : start,
#                              "end" : stop,
#                               }
#               }
# class_data['regions'].append(disordered_region)


classifier = MrClassifier()
classifier.get_predictions(seqin, topcons_dir="/opt/MrParse/data/Q13586/topcons", jpred_dir="/opt/MrParse/data/Q13586/jpred")
pfam_data = classifier.pfam_data()
pfam_data['regions'] = region_data

rdata = 'console.log("LOADED DATA");\nvar pfam_json = %s;\n' % json.dumps(pfam_data)
with open(os.path.join(html_dir, 'data.js'), 'w') as w:
    w.write(rdata)
