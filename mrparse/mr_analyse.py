'''
Created on 18 Oct 2018

@author: jmht
'''

from mrparse.mr_hkl import HklInfo
from mrparse.mr_search_model import SearchModelFinder

class DataDisplay():
    def __init__(self, hkl_info=None, search_model_finder=None):
        self.hkl_info = hkl_info
        self.search_model_finder = search_model_finder
    
    def as_html(self):        
        html = "<html>\n<body>\n"
        if self.hkl_info:
            html += "<h2>HKLIN info</h2>{}\n".format(self.hkl_info.as_dataframe().to_html())
        if self.search_model_finder:
            html += "<h2>Search Model info</h2>{}\n".format(self.search_model_finder.as_dataframe().to_html())
        html += "</html>\n</body>\n\n"
        return html
        

def run(hklin, seqin):
    hkl_info  = HklInfo(hklin)
     
    search_model_finder = SearchModelFinder(seqin)
     
    data_display = DataDisplay(hkl_info=hkl_info, search_model_finder=search_model_finder)
     
    with open('display.html', 'w') as w:
        w.write(data_display.as_html())
    return