'''
Created on 14 Nov 2018

@author: jmht
'''
import os

def parse_jpred_output(jpred_rundir):
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


def secondary_structure_prediction(jpred_rundir):
    jpred_rundir = "/opt/MrParse/data/Q13586/jpred"
    
    ss_pred, cc_28 = parse_jpred_output(jpred_rundir)
    seqlen = len(ss_pred)
    
#     print(ss_pred)
#     print(cc_28)
    
#     cc_28 = get_start_stop(cc_28.lower(), 'c')
#     cc_28 = annotate_chunks(cc_28, 'CC')
#     cc_28_d = pfam_dict(cc_28, seqlen)
#     
#     ss_pred = get_start_stop(ss_pred.lower(), 'h')
#     ss_pred = annotate_chunks(ss_pred, 'helix')
#     ss_pred_d = pfam_dict(ss_pred, seqlen)
    # print(cc_28_d)
    # print(ss_pred_d)
    
    return ss_pred


