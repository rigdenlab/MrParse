import pickle
import set_mrparse_path

from mrparse.mr_search_model import RegionFinder, find_hits, get_homologs, calculate_ellg, ellg_data_from_phaser_log

dfinder = RegionFinder()
seqin = '../data/2uvoA.fasta'
hits = find_hits(seqin)
with open('hits.pkl', 'w') as w:
    pickle.dump(hits, w)
domains = dfinder.find_domains_from_hits(hits)
with open('domains.pkl', 'w') as w:
    pickle.dump(domains, w)
homologs = get_homologs(hits, domains)
# with open('homologs.pkl', 'w') as w:
#     pickle.dump(homologs, w)
 
# with open('homologs1.pkl') as f:
#     homologs = pickle.load(f) 
# ellg_data_from_phaser_log('phaser1.log', homologs)
 
calculate_ellg(homologs)
with open('homologs.pkl', 'w') as w:
    pickle.dump(homologs, w)
        
# with open('homologs.pkl') as f:
#     homologs = pickle.load(f)
# for h, d in homologs.items():
#     print h, d

# ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py