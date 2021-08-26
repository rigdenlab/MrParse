"""
Created on 16 Nov 2018

@author: jmht
"""

from mrparse.mr_topcons import TM
from mrparse.mr_deepcoil import CC
from mrparse.mr_jpred import HELIX, SHEET
from mrparse.mr_annotation import get_annotation_chunks

import colorsys


def get_N_HexCol(N=5):
    HSV_tuples = [(x * 1.0 / N, 0.5, 0.5) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x * 255), colorsys.hsv_to_rgb(*rgb))
        hex_out.append("".join(map(lambda x: chr(x).encode('hex'), rgb)))
    return ['#' + h for h in hex_out]


def add_pfam_dict_to_homologs(homologs, sequence_length):
    # Need a better way of getting the number of regions
    nregions = len(set([h.region_id for h in homologs.values()]))
    region_colors = get_N_HexCol(nregions)
    for h in homologs.values():
        start = h.query_start
        stop = h.query_stop
        name = h.name
        region_id = h.region_id
        d = {'startStyle': "curved",
             'endStyle': "curved",
             'start': start,
             'end': stop,
             'aliStart': start,
             'aliEnd': stop,
             'colour': region_colors[h.region_index],
             'text': name,
             'metadata': {"description": "Homolog {} from region #{}".format(name, region_id),
                          "database": "PHMMER search",
                          "start": start,
                          "end": stop}
             }
        jdict = {'length': sequence_length,
                 'regions': [d]}
        h._pfam_json = jdict


def add_pfam_dict_to_models(models, sequence_length):
    # Need a better way of getting the number of regions
    for m in models.values():
        start = m.query_start
        stop = m.query_stop + 1
        name = m.name
        region_id = m.region_id
        plddt_regions = m.plddt_regions
        # 'colour': '#193f90'
        d = {'startStyle': "curved",
             'endStyle': "curved",
             'start': start,
             'end': stop,
             'aliStart': start,
             'aliEnd': stop,
             'colour': '#d3d3d3',
             'text': "",
             'metadata': {"description": "Model {} from region #{}".format(name, region_id),
                          "database": "EBI AlphaFold database",
                          "start": start,
                          "end": stop}
             }
        motifs = []
        colors = {'v_low': "ff7d45",
                  "low": "ffdb13",
                  "confident": "65cbf3",
                  "v_high": "0053d6"}
        for quality in plddt_regions.keys():
            for plddt_region in plddt_regions[quality]:
                plddt_region_start = plddt_region[0]
                plddt_region_end = plddt_region[1]
                md = {'colour': colors[quality],
                      'metadata': {"description": quality,
                                   "database": "EBI AlphaFold database",
                                   "type": "Quality",
                                   "end": plddt_region_end,
                                   "start": plddt_region_start},
                      "type": "pLDDT",
                      "display": "true",
                      "end": plddt_region_end,
                      "start": plddt_region_start
                      }
                motifs.append(md)

        jdict = {'length': sequence_length,
                 'regions': [d],
                 'motifs': motifs
                 }
        m._pfam_json = jdict


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
        elif chunk.annotation == HELIX:
            colour = "#ff0000"
            text = HELIX.name
            meta_desc = "Helix region #%d" % idx
        elif chunk.annotation == SHEET:
            colour = "#0000ff"
            text = SHEET.name
            meta_desc = "Beta strand region #%d" % idx
        d = {'startStyle': "straight",
             'endStyle': "straight",
             'start': chunk.start,
             'end': chunk.end,
             'aliStart': chunk.start,
             'aliEnd': chunk.end,
             'colour': colour,
             'text': text,
             'metadata': {"description": meta_desc,
                          "database": chunk.annotation.source,
                          "start": chunk.start,
                          "end": chunk.end,
                          }
             }
        regions.append(d)
    vis_data = {'length': seqlen,
                'regions': regions}
    return vis_data
