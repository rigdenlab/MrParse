'''
Created on 16 Nov 2018

@author: jmht
'''

from mrparse.mr_topcons import TM
from mrparse.mr_deepcoil import CC
from mrparse.mr_jpred import HELIX, SHEET
from mrparse.mr_annotation import get_annotation_chunks

import colorsys

def get_N_HexCol(N=5):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append("".join(map(lambda x: chr(x).encode('hex'), rgb)))
    return ['#' + h for h in hex_out]


def add_pfam_dict_to_homologs(regions, seqlen):
    if not (regions and len(regions)):
        raise RuntimeError("Need regions argument to contain the search model regions.")
    region_colors = get_N_HexCol(len(regions))
    for region in regions:
        for hit in region.matches:
            homolog = hit._homolog
            start = hit.query_start
            stop = hit.query_stop
            name = homolog.name
            d = { 'startStyle': "curved",
                  'endStyle': "curved",
                  'start': start,
                  'end': stop,
                  'aliStart': start,
                  'aliEnd': stop,
                  'colour': region_colors[region.index],
                  'text': name,
                  'metadata' : { "description" : "Homolog {} from region #{}".format(name, region.id),
                                 "database" : "PHMMER search",
                                 "start" : start,
                                 "end" : stop,
                                  }
                  }            
            jdict = {'length' : seqlen,
                     'regions' : [d]}
            homolog._pfam_json = jdict
    return


def pfam_dict_from_annotation(annotation):
    annotation_chunks = get_annotation_chunks(annotation)
    return pfam_dict_from_chunks(annotation_chunks, len(annotation))


def pfam_dict_from_chunks(chunk_data, seqlen):
    if not chunk_data:
        return None
    regions = []
    for i, chunk in enumerate(chunk_data):
        idx = i + 1
        if chunk.annotation == CC:
            colour = "#00ff00"
            text = CC.name
            meta_desc = "Coiled-coil region #%d" % idx
        elif chunk.annotation == TM:
            colour = "#aaaaaa"
            text = TM.name
            meta_desc = "Transmembrane region #%d" % idx
        elif chunk.annotation ==  HELIX:
            colour = "#ff0000"
            text = HELIX.name
            meta_desc = "Helix region #%d" % idx
        elif chunk.annotation == SHEET:
            colour = "#0000ff"
            text = SHEET.name
            meta_desc = "Beta strand region #%d" % idx
        d = { 'startStyle': "straight",
              'endStyle': "straight",
              'start': chunk.start,
              'end': chunk.end,
              'aliStart': chunk.start,
              'aliEnd': chunk.end,
              'colour': colour,
              'text': text,
              'metadata' : { "description" : meta_desc,
                             "database" : chunk.annotation.source,
                             "start" : chunk.start,
                             "end" : chunk.end,
                              }
              }
        regions.append(d)       
    vis_data = {'length' : seqlen,
                'regions' :regions}
    return vis_data

