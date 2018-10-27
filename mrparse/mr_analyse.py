'''
Created on 18 Oct 2018

@author: jmht
'''

from mr_hkl import HklInfo
from mr_search_model import SearchModelFinder

import json
import colorsys
def get_N_HexCol(N=5):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append("".join(map(lambda x: chr(x).encode('hex'), rgb)))
    hex_out = ['#' + h for h in hex_out]
    return hex_out

class RegionsDisplay():
    def __init__(self, seqlen, regions):
        self.seqlen = seqlen
        self.regions = regions
        self.region_colors = get_N_HexCol(len(regions))
        return
    
    def generate_pfam_json(self):
        region_data = [self.region_pfam_dict(i, region) for i, region in enumerate(self.regions)]
        return json.dumps(region_data)
    
    def region_pfam_dict(self, idx, region):
        start, stop = region.start_stop
        d = { 'startStyle': "curved",
              'endStyle':   "curved",
              'start':      start,
              'end':        stop,
              'aliStart':   start,
              'aliEnd':     stop,
              'colour':     self.region_colors[idx],
              'text':       str(region.ID) }
        jdict = {'length' : self.seqlen,
                 'regions' : [d]}
        return jdict

class DataDisplay():
    def __init__(self, hkl_info=None, search_model_finder=None):
        self.hkl_info = hkl_info
        self.search_model_finder = search_model_finder
    
    def as_html(self, fpath=None):        
        html = "<html>\n<body>\n"
        if self.hkl_info:
            html += "<h2>HKLIN info</h2>\n{}\n".format(self.hkl_info.as_html())
        if self.search_model_finder:
            html += "<h2>Search Model info</h2>\n{}\n".format(self.search_model_finder.as_html())
        html += "</html>\n</body>\n\n"
        if fpath:
            with open(fpath, 'w') as w:
                w.write(html)
        return html
        
def run(hklin, seqin):
    hkl_info  = HklInfo(hklin)
    search_model_finder = SearchModelFinder(seqin)
    data_display = DataDisplay(hkl_info=hkl_info, search_model_finder=search_model_finder)
    data_display.as_html('display.html')
    return