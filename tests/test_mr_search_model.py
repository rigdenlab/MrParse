import pickle
import set_mrparse_path

from mrparse.mr_search_model import DomainFinder, find_hits, get_homologs, calculate_ellg, ellg_data_from_phaser_log

# dfinder = DomainFinder()
# seqin = '../data/2uvoA.fasta'
# hits = find_hits(seqin)
# domains = dfinder.find_domains_from_hits(hits)
# homologs = get_homologs(hits, domains)
# with open('homologs1.pkl', 'w') as w:
#     pickle.dump(homologs, w)
 
# with open('homologs1.pkl') as f:
#     homologs = pickle.load(f) 
# ellg_data_from_phaser_log('phaser1.log', homologs)
 
# calculate_ellg(homologs)
# with open('homologs.pkl', 'w') as w:
#     pickle.dump(homologs, w)
    
    
with open('homologs.pkl') as f:
    homologs = pickle.load(f)
# for h, d in homologs.items():
#     print h, d
    
import pandas as pd
 
homologs = [h.__dict__ for _, h in homologs.items()]
df = pd.DataFrame(homologs, columns = ['name', 'eLLG', 'ncopies', 'molecular_weight', 'rmsd', 'seqid', 'frac_scat', 'total_frac_scat', 'total_frac_scat_known', 'pdb'])
df['eLLG'] = df['eLLG'].astype(float)
df.sort_values('eLLG', ascending=False, inplace=True)
with open('foo.html', 'w') as w:
    w.write(df.to_html())


# ccp4-src-2016-02-10/checkout/cctbx-phaser-dials-2015-12-22/phaser/phaser/CalcCCFromMRsolutions.py