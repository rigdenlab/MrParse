'''
Created on 26 Oct 2018

@author: jmht
'''
import pickle
import set_mrparse_path
from mrparse.mr_analyse import RegionsDisplay

seqin = '../data/2uvoA.fasta'
seqin = '../data/5u4pA.fasta'
seqin = '../data/5u4pB.fasta'
if True:
    from mrparse.mr_region import RegionFinder
    from mrparse.mr_hit import find_hits
    from ample.util.sequence_util import Sequence
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
rdata = rg.generate_pfam_json()
rdata = 'var JSON = %s    ;\n' % rdata
with open('data.json', 'w') as w:
    w.write(rdata)
