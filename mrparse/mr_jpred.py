'''
Created on 14 Nov 2018

@author: jmht
'''
import os


def get_jpred_prediction(jpred_rundir):
    if not os.path.isdir(jpred_rundir):
        raise RuntimeError("Cannot find directory:%s" % jpred_rundir)
    out_concise = [f for f in os.listdir(jpred_rundir) if f.endswith('.concise')][0]
    out_concise = os.path.join(jpred_rundir, out_concise)
    cc_28 = None
    ss_pred = None
    with open(out_concise) as f:
        line = f.readline()
        while line:
            prefix = 'Lupas_28:'
            if line.startswith(prefix):
                line = line.strip().replace(prefix,'')
                cc_28 = "".join(line.split(","))
            prefix = 'jnetpred:'
            if line.startswith(prefix):
                line = line.strip().replace(prefix,'')
                ss_pred = "".join(line.split(","))
            line = f.readline()
    assert cc_28 and ss_pred
    return ss_pred, cc_28


def get_start_stop(sequence, marker):
    assert isinstance(marker, str) and len(marker) == 1
    regions = []
    start = None
    for i, s in enumerate(sequence):
        if s == marker:
            if not start:
                start = i
        else:
            if start:
                regions.append((start, i))
                start = None
    if start:
        regions.append((start, i))
    return regions


def annotate_chunks(chunks, stype):
    return [(stype, start, stop) for (start, stop) in chunks]


def pfam_dict(chunk_data, seqlen,):
    regions = []
    for i, (stype, start, stop) in enumerate(chunk_data):
        idx = i + 1
        if stype == 'CC':
            colour = "#00ff00"
            text = 'CC'
            meta_desc = "Coiled-coil region #%d" % idx
            meta_db = "From Predictor"

        elif stype == 'helix':
            colour = "#ff0000"
            text = 'helix'
            meta_desc = "Helix region #%d" % idx
            meta_db = "From jpred"
            
        d = { 'startStyle': "straight",
              'endStyle': "straight",
              'start': start,
              'end': stop,
              'aliStart': start,
              'aliEnd': stop,
              'colour': colour,
              'text': text,
              'metadata' : { "description" : meta_desc,
                             "database" : meta_db,
                             "start" : start,
                             "end" : stop,
                              }
              }
        regions.append(d)       
    vis_data = {'length' : seqlen,
                'regions' :regions}
    return vis_data



def get_jpred_prediction_data():
    jpred_rundir = "/opt/MrParse/O75410_800"
    
    ss_pred, cc_28 = get_jpred_prediction(jpred_rundir)
    seqlen = len(ss_pred)
#     print(ss_pred)
#     print(cc_28)
    
    cc_28 = get_start_stop(cc_28.lower(), 'c')
    cc_28 = annotate_chunks(cc_28, 'CC')
    cc_28_d = pfam_dict(cc_28, seqlen)
    
    ss_pred = get_start_stop(ss_pred.lower(), 'h')
    ss_pred = annotate_chunks(ss_pred, 'helix')
    ss_pred_d = pfam_dict(ss_pred, seqlen)
    # print(cc_28_d)
    # print(ss_pred_d)
    
    return cc_28_d, ss_pred_d


