'''
Created on 26 Oct 2018

@author: jmht

https://omictools.com/protein-coiled-coil-domain-prediction-category

paircoil2 MIT
deepcoil **


'''
import json
import pickle
import set_mrparse_path
from ample.util.sequence_util import Sequence
from mrparse.mr_analyse import RegionsDisplay
from mrparse.mr_classify import coiled_coil_pfam_dict

seqin = '../data/5u4pA.fasta'
seqin = '../data/5u4pB.fasta'
seqin = '../data/O75410.fasta'
seqin = '../data/2uvoA.fasta'
seqin = '../data/O75410.fasta'
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
rdata = rg.generate_pfam_data()
ccdata = coiled_coil_pfam_dict(seqin)
gdata = {'classification' : ccdata,
         'regions' : rdata }

rdata = 'console.log("LOADED DATA");\nvar pfam_json = %s;\n' % json.dumps(gdata)
with open('data.js', 'w') as w:
    w.write(rdata)
