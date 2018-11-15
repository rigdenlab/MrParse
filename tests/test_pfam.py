'''
Created on 26 Oct 2018

@author: jmht

https://omictools.com/protein-coiled-coil-domain-prediction-category

paircoil2 MIT
deepcoil **


'''
import json
import pickle
import sys
import set_mrparse_path
from ample.util.sequence_util import Sequence
from mrparse.mr_analyse import RegionsDisplay
from mrparse.mr_classify import coiled_coil_pfam_dict
from mrparse.mr_jpred import get_jpred_prediction_data




seqin = '../data/5u4pA.fasta'
seqin = '../data/5u4pB.fasta'
seqin = '../data/O75410.fasta'
seqin = '../data/2uvoA.fasta'
seqin = '../data/O75410.fasta'
seqin = '../O75410_800/O75410_800.fasta'
if True:
    from mrparse.mr_region import RegionFinder
    from mrparse.mr_hit import find_hits
    rfinder = RegionFinder()
    hits = find_hits(seqin)
    with open('hits.pkl', 'w') as w:
        pickle.dump(hits, w)
    regions = rfinder.find_regions_from_hits(hits)
    with open('regions.pkl', 'w') as w:
        pickle.dump(regions, w)
else:
    with open('regions.pkl') as regions_fh:
        regions = pickle.load(regions_fh)
 
seqlen = len(Sequence(fasta=seqin).sequence())
rg = RegionsDisplay(seqlen, regions)
region_data = rg.generate_pfam_data()

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

class_data, sspred_data = get_jpred_prediction_data()

gdata = {'classification' : class_data,
         'ss_pred' : sspred_data,
         'regions' : region_data }

rdata = 'console.log("LOADED DATA");\nvar pfam_json = %s;\n' % json.dumps(gdata)
with open('data.js', 'w') as w:
    w.write(rdata)
