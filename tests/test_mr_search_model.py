import pickle
import set_mrparse_path

from mrparse.mr_search_model import RegionFinder
from mrparse.mr_homolog import get_homologs, calculate_ellg
from mrparse.mr_hit import find_hits


rfinder = RegionFinder()
seqin = '../data/2uvoA.fasta'
hits = find_hits(seqin)
with open('hits.pkl', 'w') as w:
    pickle.dump(hits, w)
regions = rfinder.find_regions_from_hits(hits)
with open('regions.pkl', 'w') as w:
    pickle.dump(regions, w)
    
# homologs = get_homologs(hits, regions)
# with open('homologs.pkl', 'w') as w:
#     pickle.dump(homologs, w)
 
# with open('homologs1.pkl') as f:
#     homologs = pickle.load(f) 
# ellg_data_from_phaser_log('phaser1.log', homologs)
 
# calculate_ellg(homologs)
# with open('homologs.pkl', 'w') as w:
#     pickle.dump(homologs, w)
        
# with open('homologs.pkl') as f:
#     homologs = pickle.load(f)
# for h, d in homologs.items():
#     print h, d

# ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py