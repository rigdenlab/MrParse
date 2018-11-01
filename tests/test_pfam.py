'''
Created on 26 Oct 2018

@author: jmht
'''
import sys
import pickle
import set_mrparse_path

from mrparse.mr_analyse import RegionsDisplay

with open('regions.pkl') as regions_fh:
    regions = pickle.load(regions_fh)
       
seqlen = 172
rg = RegionsDisplay(seqlen, regions)
rdata = rg.generate_pfam_json()
rdata = 'var JSON = %s    ;\n' % rdata
with open('data.json', 'w') as w:
    w.write(rdata)

# root_div = document,querySelector('div.dg')
# for i, region in enumerate(data):
#     div_id = 'dg_' + str(i)
#     gdiv = document.createElment('div', id=div_id, class='pfam_dg') # use class to fix any display issues
#     pg = new PfamGraphic();
#     pg.setParent(div_id)
#     pg.setSequence(region);
#     pg.render();