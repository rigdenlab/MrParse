'''
Created on 16 Nov 2018

@author: jmht
'''

from mrparse.mr_sequence import CC_SYMBOL, TM_SYMBOL, HELIX_SYMBOL, BSHEET_SYMBOL


import colorsys
def get_N_HexCol(N=5):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append("".join(map(lambda x: chr(x).encode('hex'), rgb)))
    hex_out = ['#' + h for h in hex_out]

    
def generate_region_pfam_data(regions, seqlen):
    region_colors = get_N_HexCol(len(regions))
    region_data = []
    for idx_region, region in enumerate(regions):
        for idx_range, rrange in enumerate(region.ranges):
            start, stop = map(int, rrange.split('-'))
            name = region.matches[idx_range]
            d = { 'startStyle': "curved",
                  'endStyle': "curved",
                  'start': start,
                  'end': stop,
                  'aliStart': start,
                  'aliEnd': stop,
                  'colour': region_colors[idx_region],
#                       'text': str(region.ID)}
                  'text': name,
                  'metadata' : { "description" : "Domain #%s" % region.ID,
                                 "database" : "PHMMER search",
                                 "start" : start,
                                 "end" : stop,
                                  }
                  }            
            jdict = {'length' : seqlen,
                     'regions' : [d]}
            region_data.append(jdict)
    return region_data

def pfam_classification_dict(chunk_data, seqlen):
    regions = []
    for i, chunk in enumerate(chunk_data):
        idx = i + 1
        if chunk.stype == 'CC':
            colour = "#00ff00"
            text = 'CC'
            meta_desc = "Coiled-coil region #%d" % idx
        elif chunk.stype == 'TM':
            colour = "#0000ff"
            text = 'TM'
            meta_desc = "Transmembrane region #%d" % idx
            
        d = { 'startStyle': "straight",
              'endStyle': "straight",
              'start': start,
              'end': stop,
              'aliStart': start,
              'aliEnd': stop,
              'colour': colour,
              'text': text,
              'metadata' : { "description" : meta_desc,
                             "database" : database,
                             "start" : start,
                             "end" : stop,
                              }
              }
        regions.append(d)       
    vis_data = {'length' : seqlen,
                'regions' :regions}
    return vis_data

def pfam_secondary_structure_dict(chunk_data, seqlen):
    regions = []
    for i, chunk in enumerate(chunk_data):
        if chunk.stype == 'helix':
            colour = "#ff0000"
            text = 'helix'
            meta_desc = "Helix region #%d" % idx
        d = { 'startStyle': "straight",
              'endStyle': "straight",
              'start': start,
              'end': stop,
              'aliStart': start,
              'aliEnd': stop,
              'colour': colour,
              'text': text,
              'metadata' : { "description" : meta_desc,
                             "database" : database,
                             "start" : start,
                             "end" : stop,
                              }
              }
        regions.append(d)       
    vis_data = {'length' : seqlen,
                'regions' :regions}
    return vis_data